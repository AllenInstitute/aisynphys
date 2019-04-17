from __future__ import print_function
import os, re, json
from . import config


_lims_engine = None
def lims_engine():
    global _lims_engine
    if _lims_engine is None:
        from sqlalchemy import create_engine
        from sqlalchemy.pool import NullPool
        _lims_engine = create_engine(config.lims_address, poolclass=NullPool)
    return _lims_engine


def query(query_str):
    conn = lims_engine().connect()
    try:
        result = conn.execute(query_str).fetchall()
    finally:
        conn.close()
    return result



def specimen_info(specimen_name=None, specimen_id=None):
    """Return a dictionary of information about a slice specimen queried from LIMS.
    
    Also generates information about the hemisphere and which side of the slice
    was patched.
    
    Returns
    -------
    organism : "mouse" or "human"
    age : age of specimen in days at time of sectioning
    date_of_birth : donor's date of birth
    genotype : the full genotype of the donor as recorded in labtracks
    weight : weight in grams
    sex : 'M', 'F', or 
    plane_of_section : 'coronal' or 'sagittal'
    hemisphere : 'left' or 'right'
    thickness : speimen slice thickness (unscaled meters)
    section_instructions : description of the slice angle and target region
        used for sectioning
    flipped : boolean; if True, then the slice was flipped relative to its 
        blockface image during recording
    sectioning_mount_side : the side of the tissue that was mounted during
        sectioning
    exposed_surface : The surface that was exposed during the experiment (right, 
        left, anterior, or posterior)
    section_number : indicates the order this slice was sectioned (1=first)
    """
    
    # Query all interesting information about this specimen from LIMS
    q = """
        select 
            organisms.name as organism, 
            ages.days as age,
            donors.date_of_birth as date_of_birth,
            donors.full_genotype as genotype,
            donors.weight as weight,
            genders.name as sex,
            tissue_processings.section_thickness_um as thickness,
            tissue_processings.instructions as section_instructions,
            plane_of_sections.name as plane_of_section,
            flipped_specimens.name as flipped,
            specimens.histology_well_name as histology_well_name,
            specimens.carousel_well_name as carousel_well_name,
            specimens.parent_id as parent_id,
            specimens.name as specimen_name,
            specimens.id as specimen_id
        from specimens
            left join donors on specimens.donor_id=donors.id 
            left join organisms on donors.organism_id=organisms.id
            left join ages on donors.age_id=ages.id
            left join genders on donors.gender_id=genders.id
            left join tissue_processings on specimens.tissue_processing_id=tissue_processings.id
            left join plane_of_sections on tissue_processings.plane_of_section_id=plane_of_sections.id
            left join flipped_specimens on flipped_specimens.id = specimens.flipped_specimen_id
    """
    if specimen_name is not None:
        sid = specimen_name.strip()
        q += "where specimens.name='%s';" % sid
    elif specimen_id is not None:
        sid = specimen_id
        q += "where specimens.id='%d';" % sid
    else:
        raise ValueError("Must specify specimen name or ID")
        
    r = query(q)
    if len(r) != 1:
        raise Exception("LIMS lookup for specimen '%s' returned %d results (expected 1)" % (sid, len(r)))
    rec = dict(r[0])
    
    # convert thickness to unscaled
    rec['thickness'] = rec['thickness'] * 1e-6
    # convert organism to more easily searchable form
    rec['organism'] = {'Mus musculus': 'mouse', 'Homo Sapiens': 'human'}[rec['organism']]
    # convert flipped to bool
    rec['flipped'] = {'flipped': True, 'not flipped': False, 'not checked': None, 'unknown': 'unknown'}.get(rec['flipped'])
    
    # Parse the specimen name to extract more information about the plane of section.
    # Mouse format is:  
    #    driver1;reporter1;driver2;reporter2-AAAAAA.BB.CC
    #        AAAAAA = donor ID
    #            BB = slice number
    #            CC = orientation and hemisphere
    spec_name = rec['specimen_name']
    if rec['organism'] == 'mouse':
        m = re.match(r'(.*)(-(\d{6,7}))?(\.(\d{2}))(\.(\d{2}))$', spec_name)
        if m is None:
            raise Exception('Could not parse mouse specimen name: "%s"' % spec_name)
        
        rec['section_number'] = int(m.groups()[4])
        
        # The last number contains information about the orientation and hemisphere
        orientation_num = m.groups()[6]
        plane, hem, mount = {
            '01': ('coronal', 'left', 'anterior'),
            '02': ('coronal', 'right', 'anterior'),
            '03': ('coronal', 'left', 'posterior'),
            '04': ('coronal', 'right', 'posterior'),
            '05': ('sagittal', 'left', 'right'),
            '06': ('sagittal', 'right', 'left'),
        }[orientation_num]
        
        if plane != rec['plane_of_section']:
            raise ValueError("Slice orientation from specimen name (.%s=%s) does not match plane_of_section recorded in LIMS (%s)" % (orientation_num, plane, rec['plane_of_section']))
        
        rec['hemisphere'] = hem
        rec['sectioning_mount_side'] = mount
        
        # decide which surface was patched
        if rec['flipped'] is True:
            rec['exposed_surface'] = mount
        elif rec['flipped'] is False:
            rec['exposed_surface'] = {'anterior': 'posterior', 'posterior': 'anterior', 'right': 'left', 'left': 'right'}[mount]
        else:    
            rec['exposed_surface'] = None
            
    # Human format is:
    #   Haa.bb.ccc.dd.ee.ff
    elif rec['organism'] == 'human':
        m = re.match(r'H(\d+)\.(\d+)\.(\d+)\.(\d+)\.(\d+)(\.(\d+))?$', spec_name)
        if m is None:
            raise Exception('Could not parse human specimen name: "%s"' % spec_name)
        rec['hemisphere'] = None
        rec['sectioning_mount_side'] = None
        rec['exposed_surface'] = None
        rec['human_donor_site'] = int(m.groups()[1])
        rec['human_donor_number'] = int(m.groups()[2])
        rec['block_number'] = int(m.groups()[3])
        rec['section_number'] = int(m.groups()[4])
        rec['subsection_number'] = None if m.groups()[6] is None else int(m.groups()[6])
        
        
    else:
        raise Exception('Unsupported organism: "%s"' % rec['organism'])
    
    return rec
    
    
def specimen_images(specimen):
    """Return a list of dicts describing images for a specimen.

    Each dict looks like::

        {'id': sub_image_id, 'file': image_file_path, 'treatment': treatment_name, 'resolution': um_per_px, 'url': lims_url}

    If an image is a stack, then the values for 'id' and 'file' will be lists.
    
    Parameters
    ----------
    specimen : int | str
        Either the ID (int) or name (str) of the specimen.
    """

    field = 'id' if isinstance(specimen, int) else 'name'
    images = []

    # First get a list of all image series for the specimen
    q = """
        select image_series.id, image_series.is_stack from specimens 
        join image_series on image_series.specimen_id=specimens.id 
        where specimens.%s='%s' and
        image_series.type='FocalPlaneImageSeries';
        """ % (field, specimen)

    # for each image series, get all sub images and decide whether to treat them as 
    # a stack or a set of images with different treatments
    for image_series in query(q):
        q = """
            select distinct sub_images.id, images.jp2, scans.resolution, treatments.name, slides.storage_directory from image_series
            join sub_images on sub_images.image_series_id=image_series.id
            join images on images.id = sub_images.image_id
            left join treatments on treatments.id = images.treatment_id
            left join slides on slides.id=images.slide_id
            left join scans on scans.slide_id=slides.id
            where image_series.id=%d;
            """ % image_series['id']
        results = query(q)

        if image_series['is_stack'] is True:
            image_ids = {}
            image_files = {}
            # sift through stack and group images by treatment and resolution
            for image in results:
                key = (image['name'], image['resolution'])
                image_ids.setdefault(key, []).append(image['id'])
                image_files.setdefault(key, []).append(image['storage_directory'].rstrip('/') + '/' + image['jp2'])
            for k in image_ids:
                # not sure how to generate an image stack url
                images.append({'id':image_ids[k], 'file': image_files[k], 'treatment': k[0], 'resolution': k[1], 'url': None, 'image_series': image_series['id']})
        else:
            for image in results:
                if image['storage_directory'] is None:
                    path = None
                else:
                    path = image['storage_directory'].rstrip('/') + '/' + image['jp2']
                url = "http://lims2/siv?sub_image=%d" % image['id']
                images.append({'id':image['id'], 'file': path, 'treatment': image['name'], 'resolution': image['resolution'], 'url': url, 'image_series': image_series['id']})
            
    return images


def specimen_20x_image(specimen, treatment='Biocytin'):
    """Return the path to a 20x aff image file.
    """
    images = specimen_images(specimen)
    for image in images:
        if image['treatment'] == treatment:
            file_base = os.path.splitext(image['file'])[0]
            # Double // ensures the path is treated as a network location on windows.
            # On posix, this should have no effect.
            return '//' + file_base.lstrip('/') + '.aff'
    return None


def specimen_species(specimen_name):
    """returns species information
    """
    q = """
    select organisms.name as species 
    from specimens left join donors on specimens.donor_id = donors.id
    left join organisms on donors.organism_id = organisms.id
    where specimens.name = '%s';
    """ % specimen_name
    r = query(q)
    if len(r) == 0:
        raise ValueError("Could not find donor information")
    return r[0]['species']


def specimen_id_from_name(spec_name):
    """Return the LIMS ID of a specimen give its name.
    """
    recs = query("select id from specimens where name='%s'" % spec_name)
    if len(recs) == 0:
        raise ValueError('No LIMS specimen named "%s"' % spec_name)
    return recs[0]['id']


def specimen_name(spec_id):
    recs = query("select name from specimens where id=%s" % spec_id)
    if len(recs) == 0:
        raise ValueError('No LIMS specimen with ID %d' % spec_id)
    return recs[0]['name']


def specimen_ephys_roi_plans(specimen):
    """Return a list of all ephys roi plans for this specimen.

    Parameters
    ----------
    specimen : int | str
        Either the ID (int) or name (str) of the specimen.
    """
    if not isinstance(specimen, int):
        specimen = specimen_id_from_name(specimen)
    
    recs = query("""
        select 
            ephys_roi_plans.id as ephys_roi_plan_id,
            ephys_specimen_roi_plans.id as ephys_specimen_roi_plan_id,
            ephys_roi_plans.name as name
        from 
            ephys_specimen_roi_plans
            join ephys_roi_plans on ephys_specimen_roi_plans.ephys_roi_plan_id=ephys_roi_plans.id
        where 
            ephys_specimen_roi_plans.specimen_id=%d
    """ % specimen)
    return recs


def cell_cluster_ids(specimen):
    """Return the IDs of all cell-cluster children of *specimen*.

    Parameters
    ----------
    specimen : int | str
        Either the ID (int) or name (str) of the specimen.
    """
    if not isinstance(specimen, int):
        specimen = specimen_id_from_name(specimen)
    q = """
    select specimens.id from specimens 
    join specimen_types_specimens on specimen_types_specimens.specimen_id=specimens.id
    join specimen_types on specimen_types.id=specimen_types_specimens.specimen_type_id
    where specimens.parent_id=%d
    and specimen_types.name='CellCluster'
    """ % specimen
    recs = query(q)
    return [rec['id'] for rec in recs]


def child_specimens(specimen):
    if not isinstance(specimen, int):
        specimen = specimen_id_from_name(specimen)
    q = """
    select id from specimens 
    where specimens.parent_id=%d
    """ % specimen
    recs = query(q)
    return [rec['id'] for rec in recs]    


def parent_specimen(specimen):
    if isinstance(specimen, int):
        q = """select parent_id from specimens where specimens.id=%d""" % specimen
    else:
        q = """select parent_id from specimens where specimens.name='%s'""" % specimen
    recs = query(q)
    return recs[0]['parent_id']


def cell_cluster_data_paths(specimen):
    if not isinstance(specimen, int):
        specimen = specimen_id_from_name(specimen)
    recs = query("""
        select ephys_roi_results.storage_directory 
        from specimens 
        join ephys_roi_results on ephys_roi_results.id=specimens.ephys_roi_result_id
        where specimens.id=%d
    """ % specimen)
    return [r['storage_directory'] for r in recs]

  
def specimen_metadata(specimen):
    if not isinstance(specimen, int):
        specimen = specimen_id_from_name(specimen)
    recs = query("select data from specimen_metadata where specimen_id=%d" % specimen)
    if len(recs) == 0:
        return None
    meta = recs[0]['data']
    if meta == '':
        return None

    if isinstance(meta, basestring):
        meta = json.loads(meta)  # unserialization corrects for a LIMS bug; we can remove this later.
    return meta

def specimen_tags(specimen):
    if not isinstance(specimen, int):
        specimen = specimen_id_from_name(specimen)
    q = """
    select name from specimen_tags
    join specimen_tags_specimens on specimen_tags_specimens.specimen_tag_id=specimen_tags.id
    where specimen_tags_specimens.specimen_id=%d""" % specimen
    recs = query(q)
    return [rec['name'] for rec in recs]


def specimen_type(specimen):
    if not isinstance(specimen, int):
        specimen = specimen_id_from_name(specimen)
    q = """
    select specimen_types.name from specimens 
    left join specimen_types_specimens on specimen_types_specimens.specimen_id=specimens.id
    left join specimen_types on specimen_types.id=specimen_types_specimens.specimen_type_id
    where specimens.id=%d
    """ % specimen
    recs = query(q)
    return recs[0]['name']


def slice_parent(spec_id):
    while True:
        typ = specimen_type(spec_id)
        if typ is None:
            return spec_id
        spec_id = parent_specimen(spec_id)


def filename_base(specimen_id, acq_timestamp):
    """Return a base filename string to be used for all LIMS uploads.
    """
    return "synphys-%d-%s" % (specimen_id, acq_timestamp)


def get_incoming_dir(specimen_name):
    """Returns the path for incoming files for each project
    """
    q = """
    select projects.incoming_directory from projects 
    left join specimens on projects.id = specimens.project_id
    where specimens.name='%s'
    """ % specimen_name
    recs = query(q)
    if recs[0]['incoming_directory'] == None:
        trigger_dir = get_trigger_dir(specimen_name)
        return os.path.dirname(os.path.dirname(trigger_dir))
    else:
        return recs[0]['incoming_directory']


def get_trigger_dir(specimen_name):
    """Returns the path for trigger files for each project
    """
    q = """
    select projects.trigger_dir from projects 
    left join specimens on projects.id = specimens.project_id
    where specimens.name  = '%s'
    """ % specimen_name
    recs = query(q)
    return recs[0]['trigger_dir']


def submit_expt(spec_name, acq_timestamp, nwb_file, json_file):
    import limstk.LIMStk as limstk
    #limstk.init_log()

    """
    This is an example of creating a trigger file and copying files for multi-patch.
    lims_scheduler_d.py should be running to receive requests from this library.

    When creating the session, HWBIgor and metadata are values you would be passing in (per confluence spec) rather than
    retrieving from lims2/ but 'id' will be the value you are passing to lims2/ to get the trigger directory and roi plans
    back from.

    These are defined in the lims_config.yml
    """
    spec_id = specimen_id_from_name(spec_name)

    # notice the "' notation when adding the filename - this is to allow the ' to appear in the trigger file
    # 'multipatch' is a lookup key for your specification but,
    # id, NWBIgor and metadata are key-value pairs you want to see in the trigger file but that do not come from lims
    
    filebase = filename_base(spec_id, acq_timestamp)
    
    incoming_files = dict(
        NWBIgor="/allen/programs/celltypes/production/incoming/mousecelltypes/%s.nwb" % filebase,
        metadata="/allen/programs/celltypes/production/incoming/mousecelltypes/%s.json" % filebase,
    )
    
    # they are also defined in the limstk_config.
    lims_session = limstk.Session(
        'multipatch',  # defined in the limstk_config file
        id=spec_id,
        NWBIgor="'" + incoming_files['NWBIgor'] + "'",
        metadata="'" + incoming_files['metadata'] + "'"
    )

    # Because there could be multiple plans returned on this query, the user has to determine which plan id is correct
    # For these situations, there are 'manual' requests like specimens_by_id and specimens_by_well_name
    resp = lims_session.request('specimens_by_id', id=spec_id)

    # This data can be manually added to the trigger data
    lims_session.trigger_data['id'] = resp['ephys_specimen_roi_plans'][0]['id']

    # enumerate the files you'd like to copy over
    lims_session.add_to_manifest(json_file, dst_filename='%s.json' % filebase)
    lims_session.add_to_manifest(nwb_file, dst_filename='%s.nwb' % filebase)

    # to finish out, schedule the session
    lims_session.commit_manifest(trigger_file='%s.mp' % filebase)

    return incoming_files


def expt_submissions(specimen, acq_timestamp):
    """Return information about the status of each submission found for an
    experiment, identified by its specimen ID and experiment uid.
    """
    if not isinstance(specimen, int):
        specimen = specimen_id_from_name(specimen)
    submissions = []
    filebase = filename_base(specimen, acq_timestamp)
    
    # Do we have incoming files?
    incoming_path = '/allen/programs/celltypes/production/incoming/mousecelltypes'
    incoming_nwb = os.path.join(incoming_path, filebase + '.nwb')
    incoming_trigger = os.path.join(incoming_path, 'trigger', '%s.mp' % filebase)
    failed_trigger = os.path.join(incoming_path, 'failed_trigger', '%s.mp' % filebase)
    
    if os.path.exists(incoming_nwb):
        if os.path.exists(incoming_trigger):
            submissions.append(("trigger pending", incoming_trigger))
        if os.path.exists(failed_trigger):
            error = open(failed_trigger+'.err', 'r').read()
            submissions.append(("trigger failed", failed_trigger, error))
            
    # Anything in LIMS already?
    cluster_ids = expt_cluster_ids(specimen, acq_timestamp)
    for cid in cluster_ids:
        data_path = cell_cluster_data_paths(cid)
        submissions.append(("succeeded", cid, data_path))
    
    return submissions
    

def expt_cluster_ids(specimen, acq_timestamp):
    """Return a list of CellCluster IDs associated with an experiment
    """
    if not isinstance(specimen, int):
        specimen = specimen_id_from_name(specimen)
    cluster_ids = cell_cluster_ids(specimen)
    cids = []
    for cid in cluster_ids:
        meta = specimen_metadata(cid)
        if meta is not None and meta['acq_timestamp'] == acq_timestamp:
            cids.append(cid)
    return cids


def cluster_cells(cluster):
    """Return information about a CellCluster's child cells.
    """
    if not isinstance(cluster, int):
        cluster = specimen_id_from_name(cluster)
    
    q = """select child.id, child.name, child.x_coord, child.y_coord, child.external_specimen_name, child.ephys_qc_result, biospecimen_polygons.polygon_id
    from specimens parent 
    left join specimens child on child.parent_id=parent.id
    left join biospecimen_polygons on biospecimen_polygons.biospecimen_id=child.id
    where parent.id=%d
    """ % cluster

    recs = query(q)
    return recs


if __name__ == '__main__':
    # testing specimen
    spec_name = "Ntsr1-Cre_GN220;Ai14-349905.03.06"
    spec_id = specimen_id_from_name(spec_name)
    cluster_ids = cell_cluster_ids(spec_id)
    print("Slice:", spec_id)
    print(specimen_info(spec_name))
    for cid in cluster_ids:
        print("  %d %s" % (cid, specimen_metadata(cid)))
    #print(cell_cluster_data_paths(cluster_ids[0]))
    print("")
    for sub in expt_submissions(spec_id, 1505768693.087):
        print(sub)

    #submit_expt(spec_id, 'lims_test.nwb', 'lims_test.json')
    
    
