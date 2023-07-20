import os, yaml, re, json
from datetime import datetime, timedelta
from collections import OrderedDict
import acq4

from . import config, lims, constants, genotypes
from . import yaml_local  # adds support for OrderedDict


class ExperimentMetadataSubmission(object):
    """Handles final metadata QC before submission.
    
    Outputs:
    
    - pipettes.yml - list of pipettes with patch status
    - file_manifest.yml - list of files selected for LIMS submission
    - acq4 .index files updated with file categories, LIMS specimen data, rig name.

    If possible, some fields in pipettes.yml are filled from old-format metadata files
    """
    message = "Copying submission metadata back to index files"

    _expt_list = None  # used for importing old metadata
    
    def __init__(self, site_dh, files, pipettes):
        self.site_dh = site_dh
        self.files = files
        self.pipettes = pipettes

    def check(self):
        """
        * error if any files do not exist
        * warning if existing files are not mentioned in file list
        * warning if any metadata would be overwritten with a new value
        """

        # Projects for which we do not require a LIMS specimen 
        no_lims_specimen_projects = ["human_mFISH_pilot"]

        # Projects for which we can skip dissection timing checks
        no_dissection_time_projects = ["human_mFISH_pilot"]

        errors = []
        warnings = []
        
        site_info = self.site_dh.info()
        slice_dh = self.site_dh.parent()
        slice_info = slice_dh.info()
        expt_dh = slice_dh.parent()
        expt_info = expt_dh.info()
        
        spec_name = slice_info['specimen_ID'].strip()

        # slice project is specified
        project = slice_info.get('project', '')
        if project == '':
            errors.append("Must specify slice.project")

        if project in no_lims_specimen_projects:
            # Some projects may use specimens that are not tracked in LIMS
            self.spec_info = None
        else:
            self.spec_info = lims.specimen_info(spec_name)

        
        # Do all categorized files actually exist?
        # Was a file categorized differently on a previous submission?
        # Are file sizes appropriate?
        max_sizes = {
            'MIES physiology': 5000,
            'Multipatch log': 1,
            'slice anatomy': 10,
            'slice quality stack': 5000,
            'recording site': 5000,
        }
        for file_info in self.files:
            fh = slice_dh[file_info['path']]
            if not fh.isFile():
                errors.append("file %s does not exist" % fh.name())
                continue
            info = fh.info()
            if 'category' in info and info['category'] != file_info['category']:
                warnings.append("category for file %s will change %s => %s" % 
                                (fh.name(), info['category'], file_info['category']))
            
            max_size = max_sizes.get(file_info['category'], 0)
            size = int(os.stat(fh.name()).st_size / 1e6)
            if max_size > 0 and size > max_size:
                errors.append('File "%s" is too large (%dMB > %dMB max for category "%s")' % 
                              (fh.name(), size, max_size, file_info['category']))
            
        
        # Did we forget to categorize a file?
        all_files = [f['path'] for f in self.files]
        for dh in [slice_dh, self.site_dh]:
            for fname in dh.ls():
                fh = dh[fname]
                relname = fh.name(relativeTo=slice_dh)
                if fh.isFile() and relname not in all_files and fh.ext() not in ['.pxp']:
                    warnings.append("file %s is not mentioned in metadata." % relname)
        
        # Do we already have pipette information submitted?
        self.pip_file = self.site_dh['pipettes.yml'].name()
        if os.path.exists(self.pip_file):
            yaml_pips = yaml.load(open(self.pip_file, 'rb'))
            if yaml_pips != self.pipettes:
                warnings.append('pipette metadata file will be overwritten.')
        
        # make sure we have at most one nwb file to submit
        source_nwb_files = [f['path'] for f in self.files if f['category'] == 'MIES physiology']
        have_nwb = len(source_nwb_files) > 0
        if len(source_nwb_files) == 0:
            warnings.append("No NWB files specified")
        if len(source_nwb_files) > 1:
            errors.append("%d NWB files specified (should be 0 or 1)" % len(source_nwb_files))
        
        # make sure we have one multipatch log file to submit if an nwb is present
        # or 0-1 log files if there is no nwb
        source_log_files = [f['path'] for f in self.files if f['category'] == 'Multipatch log']
        if have_nwb and len(source_log_files) == 0:
            # Occasionally user forgets to record this data; usually we can work with it
            warnings.append("No MultiPatch log files specified")
        if len(source_log_files) > 1:
            errors.append("%d MultiPatch log files specified (should be 0 or 1)" % len(source_log_files))

        
        # Sanity checks on experiment metadata:
        
        
        # Is rig name set and known?
        if config.rig_name not in ('MP0', 'MP1', 'MP2', 'MP3', 'MP4', 'MP5'):
            errors.append('Unrecognized rig name "%s"' % config.rig_name)
        
        # Is ACSF set?
        acsf = expt_info.get('solution', None)
        if acsf not in constants.ACSF_RECIPES:
            errors.append('Unrecognized ACSF recipe: "%s"' % acsf)
        
        # Is temperature set?
        temp = expt_info.get('temperature', '')
        m = re.match(r'(rt|\d{2}\s*c?)', temp.lower())
        if m is None:
            errors.append('Unrecognized temperature: "%s"' % temp)
        
        # is operator name set?
        operator = expt_info.get('rig_operator', '')
        if operator is None or not operator.startswith('Operator '):
            errors.append('Rig operator field is not set.')
        
        # Slice time ok?
        site_date = datetime.fromtimestamp(site_info['__timestamp__'])
        tod = expt_info.get('time_of_dissection', '')
        if project in no_dissection_time_projects:
            pass
        elif tod == '':
            warnings.append("Time of dissection not specified")
        else:
            m = re.match(r'((20\d{2})-(\d{1,2})-(\d{1,2}) )?(\d{1,2}):(\d{1,2})', tod)
            if m is None:
                errors.append('Time of dissection must be "HH:MM" or "YYYY-MM-DD HH:MM"')
            else:
                _, year, mon, day, h, m = m.groups()
                if int(h) < 7:
                    warnings.append("Dissection time appears to be very early. Did you mean %d:%s?" % (int(h)+12, m))
                # interpret time of dissection
                if year is None:
                    diss_time = datetime(site_date.year, site_date.month, site_date.day, int(h), int(m))
                    extra_err = " (try specifying the full date like `yyyy-mm-dd hh:mm`)"
                else:
                    diss_time = datetime(int(year), int(mon), int(day), int(h), int(m))
                    extra_err = ""
                # check time
                seconds_since_dissection = (site_date - diss_time).total_seconds()
                if seconds_since_dissection < 30*60:
                    warnings.append("Time of dissection is later than experiment time - 30 minutes")
                if seconds_since_dissection < 0:
                    errors.append("Time of dissection is later than experiment time" + extra_err)
                if seconds_since_dissection > 10*3600:
                    warnings.append("Time of dissection is more than 10 hours prior to experiment")

        # check specimen age
        if self.spec_info is not None and self.spec_info['age'] is None:
            warnings.append("Donor age is not set in LIMS.")

        # Check specimen death was today
        if self.spec_info is not None and self.spec_info['date_of_birth'] is not None:
            dod = self.spec_info['date_of_birth'].date() + timedelta(self.spec_info['age'])
            days_since_death = (site_date.date() - dod).days
            if days_since_death > 0:
                warnings.append("Specimen was dissected before today, likely need to update LabTracks info in LIMS(Date of experiment: %s, Date of death: %s, you decide)" % (site_date.date(), dod))
            if days_since_death < 0:
                warnings.append("Specimen is from the future, likely need to update LabTracks info in LIMS")

        # Sanity checks on pipette metadata:
        patchseq_tube_ids = []
        for pip_id, pip in self.pipettes.items():
            # Check pipette status
            if pip['pipette_status'] not in ['No seal', 'Low seal', 'GOhm seal', 'Technical failure', 'No attempt', 'Not recorded']:
                warnings.append('Pipette %d has unrecognized status "%s"' % (pip_id, pip['pipette_status']))
            
            # Check that patchseq tubes have correct format
            hs_name = 'HS%d' % pip_id
            tube_name = site_info.get('headstages', {}).get(hs_name, {}).get('Tube ID', '')
            if tube_name != '':
                patchseq_tube_ids.append(tube_name)
                name_check = re.match(r'P(M|T)S4_(?P<date>\d{6})_(?P<tube_id>\d{3})_A01', tube_name)
                if name_check is None:
                    errors.append('Tube ID %s for pipette %d does not have proper format' % (tube_name, pip_id))
                if name_check.group('date') != datetime.strftime(site_date, "%y%m%d"):
                    errors.append('Incorrect date in tube ID %s for pipette %d' % (tube_name, pip_id))

            # Following checks only apply if we got data from this pipette.    
            if not pip['got_data']:
                continue
            
            # Did we specify a known dye for each pipette?
            if pip['internal_dye'] not in constants.INTERNAL_DYES:
                errors.append('Pipette %d has unrecognized dye "%s"' % (pip_id, pip['internal_dye']))
            
            # Did we specify a known internal for each pipette?
            if pip['internal_solution'] not in constants.INTERNAL_RECIPES:
                errors.append('Pipette %d has unrecognized internal "%s"' % (pip_id, pip['internal_solution']))
            
            # Does the selected dye overlap with cre reporters?

        # If there are patchseq tubes make sure all IDs are unique
        if len(patchseq_tube_ids) > 0:
            if len(patchseq_tube_ids) != len(set(patchseq_tube_ids)):
                errors.append('Looks like there are duplicate patchseq Tube IDs, please double check')
        
        if project not in no_lims_specimen_projects:
            # If slice was not fixed, don't attempt LIMS submission
            try:
                sid = lims.specimen_id_from_name(spec_name)
            except ValueError as err:
                errors.append(err.message)
                sid = None
            if slice_info['plate_well_ID'] != 'not fixed' and sid is not None:
                self._check_lims(errors, warnings, spec_name, sid, site_info, slice_info)
        
        return errors, warnings
        
    def _check_lims(self, errors, warnings, spec_name, sid, site_info, slice_info):
        # LIMS upload will fail if the specimen has not been given an ephys roi plan.
        accepted_roi_plans = {
            'mouse': 'Synaptic Physiology ROI Plan', 
            'human': 'Synaptic Physiology Human ROI Plan'
        }[self.spec_info['organism']]
        roi_plans = lims.specimen_ephys_roi_plans(spec_name)
        lims_edit_href = '<a href="http://lims2/specimens/{sid}/edit">http://lims2/specimens/{sid}/edit</a>'.format(sid=sid)

        site_date = datetime.fromtimestamp(site_info['__timestamp__'])
        if site_date >= datetime(2017, 10, 1):  # newer experiments need extra checks on LIMS structures
            if len(roi_plans) == 0:
                errors.append('Specimen has no ephys roi plan. Edit:' + lims_edit_href)
            elif len(roi_plans) == 1 and roi_plans[0]['name'] not in accepted_roi_plans:
                errors.append('Specimen has wrong ephys roi plan '
                    '(expected "Synaptic Physiology ROI Plan"). Edit:' + lims_edit_href)
            elif len(roi_plans) > 1:
                roi_plan_names = [roi['name'] for roi in roi_plans]
                if accepted_roi_plans in roi_plan_names:
                    warnings.append('Specimen has multiple roi plans but Synaptic Physiology ROI plan is one of them '
                        'this is safe to submit')
                else:
                    errors.append('Specimen has multiple ephys roi plans and Synaptic Physiology ROI plan is not '
                    'one of them. Edit:' + lims_edit_href)

            # Check LIMS specimen has flipped field set
            if self.spec_info['flipped'] not in (True, False):
                if self.spec_info['organism'] == 'human':
                    # human specimens can be symmetrical enough that "flipped" is meaningless
                    warnings.append("Flipped field was not set for this specimen: %s" % lims_edit_href)
                else:
                    errors.append("Must set flipped field for this specimen: %s" % lims_edit_href)
        
        # histology well name should look something like "multi_170911_21_A01"
        hist_well = self.spec_info['histology_well_name']
        if hist_well is None:
            errors.append("Missing histology well name in LIMS: %s" % lims_edit_href)        
        else:
            if site_date > datetime(2017, 4, 30):
               # 12-well plates
                if site_date > datetime(2020, 12, 6):
                    # added M or H identifier to plate for mouse and human
                    well_name_regex_match = re.match(r'(multi|Ephys)_(\d{2})(\d{2})(\d{2})_(2[1-9])_([A-C]0[1-4])_(M|H)', hist_well)
                else:
                    well_name_regex_match = re.match(r'(multi|Ephys)_(\d{2})(\d{2})(\d{2})_(2[1-9])_([A-C]0[1-4])', hist_well)
            else:
               # older experiments used 24-well plates
               well_name_regex_match = re.match(r'(multi|Ephys)_(\d{2})(\d{2})(\d{2})_(2[1-9])_([A-D]0[1-6])', hist_well)            
            if well_name_regex_match is None:
                warnings.append("Histology well name appears to be incorrectly formatted: %s" % lims_edit_href)        
            else:
                _, yy, mm, dd, plate_n, well = well_name_regex_match.groups()[:5]
                plate_date = datetime(2000+int(yy), int(mm), int(dd))
                # find the most recent Monday
                last_monday  = site_date - timedelta(days=site_date.weekday())
                expected_plate_date = [last_monday]
                # at times of reduced staffing the plate schedule can change, add other acceptable plate dates here
                expected_plate_date.append(last_monday + timedelta(days=7)) # "next Monday"
                if abs((expected_plate_date[0] - plate_date).total_seconds()) > 7*24*3600:
                    # error if more than a week out of sync
                    errors.append("Histology well date is %s%s%s; expected %s: %s" % (yy, mm, dd, expected_plate_date.strftime('%y%m%d'), lims_edit_href))
                expected_plate_dates = [epd.date() for epd in expected_plate_date]
                if plate_date.date() not in expected_plate_dates:
                    # warning if the date is not exactly as expected
                    warnings.append("Histology well date is %s%s%s; expected %s: %s" % (yy, mm, dd, ', '.join([epd.strftime('%y%m%d') for epd in expected_plate_date]), lims_edit_href))
                    
                if int(plate_n) > 24:
                    warnings.append("Histology plate number %s is probably too high. %s" % (plate_n, lims_edit_href))
                elif int(plate_n) > 22:
                    warnings.append("Histology plate number %s might be too high? %s" % (plate_n, lims_edit_href))

        # Check carousel ID matches the one in LIMS
        cw_name = self.spec_info['carousel_well_name']
        if cw_name is None and self.spec_info['organism'] == 'human' and self.spec_info['subsection_number'] is not None:
                # this specimen was subdivided; we need to ask about the parent carousel well name instead
                # note that this behavior has changed-- newer subdivided specimens will have their own carousel well name
                parent_spec_info = lims.specimen_info(specimen_id=self.spec_info['parent_id'])
                cw_name = parent_spec_info['carousel_well_name']
        if cw_name is None:
            errors.append('No LIMS carousel well name for this specimen!')
        else:
            if cw_name != slice_info['carousel_well_ID'].strip():
                errors.append('LIMS carousel well name "%s" does not match ACQ4 carousel_well_ID "%s"' 
                    % (cw_name, slice_info['carousel_well_ID']))
            
        # If histology well name was recorded in ACQ4, make sure it matches the one in LIMS
        acq_plate_well = slice_info.get('plate_well_ID', None)
        if acq_plate_well is not None and acq_plate_well.strip() != hist_well:
            errors.append('LIMS histology well name "%s" does not match ACQ4 plate_well_ID "%s"' 
                    % (hist_well, acq_plate_well))
                    
        # make sure genotype matches specimen name
        if self.spec_info['organism'] == 'mouse':
            if self.spec_info['genotype'] is None:
                errors.append('Specimen %s has no genotype' % spec_name)
            else:
                gt = genotypes.Genotype(self.spec_info['genotype'])
                for part in gt.driver_lines + gt.reporter_lines:
                    # Error if the first chunk (usually something like "Ai193" or "Sim1") is not in the specimen name
                    part = part.partition('(')[0].partition('-')[0]
                    if part not in spec_name:
                        errors.append('Specimen name %s does not contain genotype part %s' % (spec_name, part))

    def lims_cells(self):
        reporter = {'-': 'negative', None: 'not applicable'}
        cells = []
        site_info = self.site_dh.info()
        headstages = site_info.get('headstages')
        if headstages is None:
            return
        day_info = self.site_dh.parent().parent().info()
        region = day_info.get('target_region')
        structure = 'VISp' if region == 'V1' else None
        for hs, info in headstages.items():
            cell = {
            'external_specimen_name': hs[-1],
            'patched_cell_container': None,
            'cell_reporter': None,
            'structure': None,
            }
            tube_id = info['Tube ID']
            if tube_id == '':
                cells.append(cell)
                continue
            cell['structure'] = structure
            cell['patched_cell_container'] = tube_id
            cell_reporter = info.get('Reporter')
            if cell_reporter == '-':
                cell['cell_reporter'] = 'negative'
            elif cell_reporter in ['red', 'green', 'yellow']:
                cell['cell_reporter'] = 'positive'
            else:
                cell['cell_reporter'] = 'not applicable'

            cells.append(cell)

        return cells

    def summary(self):
        summ = OrderedDict()
        summ['file categories'] = self.files
        summ['pipettes'] = self.pipettes
        
        return {'metadata': summ}
    
    def submit(self):
        slice_dh = self.site_dh.parent()
        expt_dh = slice_dh.parent()
        
        # Set rig name from config if it has not already been written
        if 'rig_name' not in expt_dh.info():
            expt_dh.setInfo(rig_name=config.rig_name)
            
        # Write LIMS info to top-level and slice directory handles.
        # This is just for convenience when browsing the raw data.
        if self.spec_info is not None:
            slice_info = self.spec_info.copy()
            donor_keys = ['organism', 'age', 'date_of_birth', 'genotype', 'weight', 'sex']
            donor_info = {k:slice_info.pop(k) for k in donor_keys}
            expt_dh.setInfo(LIMS_donor_info=donor_info)
            slice_dh.setInfo(LIMS_specimen_info=slice_info)
            
        # Write category labels for each image
        for file_info in self.files:
            if file_info['category'] == 'ignore':
                continue
            fh = slice_dh[file_info['path']]
            fh.setInfo(category=file_info['category'])
        
        # Generate yml file with pipette/cell information. Operators can edit this
        # information later on..
        yaml.dump(self.pipettes, open(self.pip_file, 'wb'), default_flow_style=False, indent=4)

        # Write an manifest file describing the data that should be uploaded to LIMS
        manifest_file = self.site_dh['file_manifest.yml'].name()
        yaml.dump(self.files, open(manifest_file, 'wb'), default_flow_style=False, indent=4)

        # Generate json for LIMS cell specimens
        cells = self.lims_cells()
        json_file = os.path.join(self.site_dh.name(), 'lims_cells.json')
        json.dump(cells, open(json_file, 'wb'))
