
from abaqus import session, getInputs
import abaqusConstants as ac
import visualization
# from caeModules import *
import numpy as np
import math, os, sys, re

def process_outputs(odb, job_name, cutoff_freq=None):
    rf_name = job_name+'-RF'
    rf_name_bot = rf_name+'-bottom'
    rf_filt_name = rf_name+'filt'
    rf_residual_name = rf_name+'-resid'
    rf_residuanorm_name = rf_name+'-residnorm'
    u_name = job_name+'-U-top'
    u_name_bot = job_name+'-U-bottom'
    disp_name = job_name+'-displacement'
    ld_name = job_name
    session.XYDataFromHistory(odb=odb, name=rf_name, outputVariableName='Reaction force: RF3 at Node 9999998 in NSET T_NODE', steps=['Load',])
    session.XYDataFromHistory(odb=odb, name=u_name, outputVariableName='Spatial displacement: U3 at Node 9999998 in NSET T_NODE', steps=['Load',])
    displacement = session.xyDataObjects[u_name]
    load = session.xyDataObjects[rf_name]
    # Get displacement at bottom node if available
    try:
        session.XYDataFromHistory(odb=odb, name=u_name_bot, outputVariableName='Spatial displacement: U3 at Node 9999999 in NSET B_NODE', steps=['Load',])
        displacement = session.xyDataObjects[u_name] - session.xyDataObjects[u_name_bot]
        displacement.setValues(sourceDescription='Displacement T_NODE - B_NODE')
    except:
        displacement = displacement * 2
        displacement.setValues(sourceDescription='Displacement T_NODE*2')
    try:
        session.XYDataFromHistory(odb=odb, name=rf_name_bot, outputVariableName='Reaction force: RF3 at Node 9999999 in NSET B_NODE', steps=['Load',])
        force_residual = session.xyDataObjects[rf_name] - session.xyDataObjects[rf_name_bot]
        session.xyDataObjects.changeKey(force_residual.name, rf_residual_name)
        force_residual.setValues(sourceDescription='Force residual T_NODE - B_NODE')
        force_residnorm = session.xyDataObjects[rf_residual_name]/np.asarray(load)[:,1].max()
        session.xyDataObjects.changeKey(force_residnorm.name, rf_residuanorm_name)
    except:
        pass
    session.xyDataObjects.changeKey(displacement.name, disp_name)
    displacement = session.xyDataObjects[disp_name]
    duration = load[-1][0]
    if cutoff_freq:
        freq = cutoff_freq
    else:
        freq = 200/duration
    xy2 = butterworthFilter(xyData=load, cutoffFrequency=freq)
    xy2.setValues(sourceDescription='butterworthFilter(cutoffFrequency={})'.format(freq))
    session.xyDataObjects.changeKey(xy2.name, rf_filt_name)
    loadfilt = session.xyDataObjects[rf_filt_name]
    xy3 = combine(displacement, loadfilt)
    session.xyDataObjects.changeKey(xy3.name, ld_name)
    print('# of history output points: {}'.format(len(displacement)))
    try:
        xy = session.XYDataFromHistory(odb=odb, name='_temp_1', outputVariableName='Percent change in mass: DMASS for Whole Model', steps=['Load',])
        print('Initial DMASS: {}'.format(xy[0][1]))
        allke = session.XYDataFromHistory(odb=odb, name=job_name+'-KE', outputVariableName='Kinetic energy: ALLKE for Whole Model', steps=['Load',])
        print('max(ALLKE) = {}'.format(currentMax( xy)[-1][1]))
        allse = session.XYDataFromHistory(odb=odb, name=job_name+'-SE', outputVariableName='Strain energy: ALLSE for Whole Model', steps=['Load',])
        kefrac = allke[10:]/allse[10:]  # Ignore first few points because there is a spike on the start of the step
        session.xyDataObjects.changeKey(kefrac.name, job_name+'-KE-div-SE')
        # session.xyDataObjectes[job_name+'-KE-div-SE']
    except:
        pass
    

def get_crack_legnths(odb, job_name, method='path', instance_name='PART-1-1', step_name='Load', path_name='Path-1',
                    coh_path_start_nset='COH_PATH_START', coh_path_end_nset='COH_PATH_END', elset_for_field='PART-1-1.COH_CL_EXTRACT',
                    a0=50.8):
    '''
    helps to extract crack lengths from the odb
    note method='field' is not working well, use method='path'
    '''

    assert(method in ('field', 'path'))
    
    dmg_variable_names = ('SDV_COH_dmg', 'SDV_CDM_d2')
    for dmg_variable_name in dmg_variable_names:
        if dmg_variable_name in odb.steps['Load'].frames[-1].fieldOutputs.keys():
            break
    else:
        print('Neither `SDV_COH_dmg` nor `SDV_CDM_d2` available in field outputs, skipping crack length measurement')
        return

    vp = session.viewports['Viewport: 1']
    vp.setValues(displayedObject=odb)
    crack_tip_dmg0 = []
    crack_tip_dmg0_all = []
    crack_tip_dmg1 = []
    length_fpz = []

    if method == 'field':
        def update_or_append(list, new_pt_time, new_pt_coord, tol=1e-5):
            updated_existing = False
            for i, pt in enumerate(list):
                if abs(new_pt_time-pt[0]) < tol:  # new_pt_time occurs at same time as previous integration point
                    if new_pt_coord > pt[1]:  # more crack extension than previously saved point, so update
                        list[i] = (new_pt_time, new_pt_coord)
                        updated_existing = True
            if not updated_existing:
                list.append((new_pt_time, new_pt_coord))
        
        session.xyDataListFromField(odb=odb, outputPosition=ac.INTEGRATION_POINT, elementSets=(elset_for_field, ),
            variable=(('COORD', ac.INTEGRATION_POINT, ((ac.COMPONENT, 'COORD1'), )), (dmg_variable_name, ac.INTEGRATION_POINT), ), )
        for xyname in session.xyDataObjects.keys():
            if xyname.startswith(dmg_variable_name):
                location = xyname.strip(dmg_variable_name).strip()
            else:
                continue
            coord1 = session.xyDataObjects['COORD:COORD1 '+location][0][1] - a0
            xyn = np.array(session.xyDataObjects[xyname])
            if any(xyn[:,1] != 0):  # Only consider integration points that experience > 0 damage at some point in history
                last_zero = xyn[np.where(xyn[:,1] == 0)[0][-1],0]
                update_or_append(crack_tip_dmg0, last_zero, coord1)
                crack_tip_dmg0_all.append((last_zero, coord1))

            if any(xyn[:,1] >= 1):
                first_one = xyn[np.where(xyn[:,1] >= 1)[0][0],0]
                # crack_tip_dmg1.append((first_one, coord1))
                update_or_append(crack_tip_dmg1, first_one, coord1)

        # TODO: sort by coord

        # Make sure there is a value for each frame
        dmg0_pts_to_add = []
        dmg1_pts_to_add = []
        for frame in odb.steps['Load'].frames:
            if frame.frameValue < crack_tip_dmg0[0][0]:
                dmg0_pts_to_add.append((frame.frameValue, 0))
            if frame.frameValue < crack_tip_dmg1[0][0]:
                dmg1_pts_to_add.append((frame.frameValue, 0))
        crack_tip_dmg0 = dmg0_pts_to_add + crack_tip_dmg0
        crack_tip_dmg1 = dmg1_pts_to_add + crack_tip_dmg1

        # Length of the fracture process zone
        tol = 1e-4
        times = [frame.frameValue for frame in odb.steps['Load'].frames]
        for time in times:
            dmg0_at_time = [x for x in crack_tip_dmg0 if abs(time-x[0]) < tol]
            dmg1_at_time = [x for x in crack_tip_dmg1 if abs(time-x[0]) < tol]
            if len(dmg0_at_time) == 1 and len(dmg1_at_time) == 1:
                lfpz = dmg0_at_time[0][1] - dmg1_at_time[0][1]
                length_fpz.append((time, lfpz))
            if len(dmg0_at_time) > 1:
                print([x for x in crack_tip_dmg0 if abs(time-x[0]) < tol])
                raise Exception('error0')
            if len(dmg1_at_time) > 1:
                raise Exception('error1')
        
        # Remove the temporary xy data
        temp_xy = [xy for xy in session.xyDataObjects.keys() if xy.startswith(dmg_variable_name) or xy.startswith('COORD:COORD1')]
        for xy in temp_xy:
            del session.xyDataObjects[xy]

    elif method == 'path':
        instance_nsets = odb.rootAssembly.instances[instance_name].nodeSets
        if coh_path_start_nset in instance_nsets and coh_path_end_nset in instance_nsets:
            print('Creating path for crack length extraction from nsets')
            start_coord = instance_nsets[coh_path_start_nset].nodes[0].coordinates
            end_coord = instance_nsets[coh_path_end_nset].nodes[0].coordinates
        else:
            start_coord = (a0, 0.5, 3.8)
            end_coord = (120.0, 0.5, 3.8)
        session.Path(name=path_name, type=ac.POINT_LIST, expression=(start_coord, end_coord))
        
        for frame_num in range(len(odb.steps[step_name].frames)):
            vp.odbDisplay.setFrame(step=0, frame=frame_num)
            vp.odbDisplay.setPrimaryVariable(variableLabel=dmg_variable_name, outputPosition=ac.INTEGRATION_POINT, )
            xy = session.XYDataFromPath(path=session.paths['Path-1'], includeIntersections=True, name='_path-dmg-frame-{}'.format(frame_num),
                shape=ac.UNDEFORMED, pathStyle=ac.PATH_POINTS, numIntervals=10, labelType=ac.X_COORDINATE, removeDuplicateXYPairs=True)
            xy = np.array(xy)
            time = odb.steps[step_name].frames[frame_num].frameValue
            if any(xy[:,1] > 0):
                ct_dmg0_idx = np.where(xy[:,1] == 0)[0][0]
                crack_tip_dmg0.append((time, xy[ct_dmg0_idx,0] - xy[0,0]))
            else:
                crack_tip_dmg0.append((time, 0))
            if any(xy[:,1] == 1.0):
                ct_dmg1_idx = np.where(xy[:,1] == 1)[0][-1]
                crack_tip_dmg1.append((time, xy[ct_dmg1_idx,0] - xy[0,0]))
                length_fpz.append((time, crack_tip_dmg0[-1][1] - crack_tip_dmg1[-1][1]))
            else:
                crack_tip_dmg1.append((time, 0))
                length_fpz.append((time, 0))


        # Remove the temporary xy data
        temp_xy = [xy for xy in session.xyDataObjects.keys() if xy.startswith('_')]
        for xy in temp_xy:
            del session.xyDataObjects[xy]
    
    session.XYData(name=job_name+'-da-dmg=0', data=crack_tip_dmg0, axis1QuantityType=visualization.QuantityType(type=ac.TIME),
                    axis2QuantityType=visualization.QuantityType(type=ac.LENGTH))
    session.XYData(name=job_name+'-da-dmg=1', data=crack_tip_dmg1, axis1QuantityType=visualization.QuantityType(type=ac.TIME),
                    axis2QuantityType=visualization.QuantityType(type=ac.LENGTH))
    session.XYData(name=job_name+'-length-fpz', data=length_fpz, axis1QuantityType=visualization.QuantityType(type=ac.TIME),
                    axis2QuantityType=visualization.QuantityType(type=ac.LENGTH))


def get_parameters(parameters_file_path):
    with open(parameters_file_path, 'r') as f:
        lines = f.readlines()

    # Definition of variables to read from paraemters file
    
    variables = [
        {'name': 'E11', 'regex': re.compile(r'^E11 .*= ([0-9.e]+)')},
        {'name': 'E22', 'regex': re.compile(r'^E22 .*= ([0-9.e]+)')},
        {'name': 'G13', 'regex': re.compile(r'^G12 .*= ([0-9.e]+)')},
        {'name': 'GIc', 'regex': re.compile(r'^GIc .*= ([0-9.e]+)')},
        {'name': 'a0','regex': re.compile(r'^crackLength .*= ([0-9.e]+)')},
        {'name': 'b','regex': re.compile(r'^width .*= ([0-9.e]+)')},
        {'name': 'h','regex': re.compile(r'^TTL .*= ([0-9.e]+)')},
    ]
    # Parse the parameters file
    param_dict = {}
    for line in lines:
        matches = [var['regex'].search(line) for var in variables]
        if any(matches):
            for i, match in enumerate(matches):
                if match:
                    name = variables[i]['name']
                    # if len(line.split('=')) > 2:
                    #     print(name, line)
                    #     print(match.group(1))
                    # else:
                    param_dict[name] = float(match.group(1))
    return param_dict

def mbt(h, b, a0, af, E11, E22, G13, GIc, npts=50, name='MBT', linear_segment='None', job_name=None):
    '''
    when linear_segment = 'from_analysis', job_name must be specified
    '''
    a = np.linspace(0.9*a0, af, num=npts)
    gamma = 1.18*np.sqrt(E11*E22)/G13
    chi = np.sqrt(E11/(11*G13)*(3-2*(gamma/(1+gamma)**2)))
    # print(chi)
    I = b*h**3/12
    force = np.sqrt((GIc*b*E11*I)/(a+chi*h)**2)
    disp = 2*force*(a+chi*h)**3/(3*E11*I)
    if linear_segment.lower() == 'connected':
        xy = np.vstack(([[0,0],], np.stack((disp, force)).T)).tolist()
    elif linear_segment.lower() == 'none':
        xy = np.stack((disp, force)).T.tolist()
    elif linear_segment.lower() == 'from_analysis':
        xy = np.stack((disp, force)).T.tolist()
        # Add separate line for initial stiffness
        # Get load-displacement from job
        ld = session.xyDataObjects[job_name]
        ld_decimated = np.array([ld[x] for x in np.linspace(0,len(ld),50).astype(int)[:-1]])
        lin_fit_data = np.array(ld_decimated[0:5])
        lin_fit_coeff = np.polyfit(lin_fit_data[:,0], lin_fit_data[:,1], 1)
        p = np.poly1d(lin_fit_coeff)
        lin_fit_x = np.array([0, ld_decimated[ld_decimated[:,1].argmax(),0]*1.1])
        print(lin_fit_x)
        lin_data = np.vstack((lin_fit_x, p(lin_fit_x))).T.tolist()
        session.XYData(data=lin_data, name=name+'_linear', axis1QuantityType=visualization.QuantityType(type=ac.DISPLACEMENT),
                    axis2QuantityType=visualization.QuantityType(type=ac.FORCE)) 
    else:
        raise Exception('Unrecognized value {} for argument linear_segement, expected one of {}'.format(
            linear_segment, ('connected', 'none', 'from_analysis')))
    session.XYData(data=xy, name=name, axis1QuantityType=visualization.QuantityType(type=ac.DISPLACEMENT),
                 axis2QuantityType=visualization.QuantityType(type=ac.FORCE))

def find_common_string(list_of_strings):
    if not list_of_strings:
        return ""
    shortest_string = min(list_of_strings, key=len)
    for i in range(len(shortest_string), -1, -1):
        for j in range(len(shortest_string) - i + 1):
            substring = shortest_string[j:j+i]
            if all(substring in s for s in list_of_strings):
                return substring
    return ""

if __name__ == '__main__':

    cutoff_freq = None  # None = default 200/step_time

    odb_names = {os.path.splitext(os.path.basename(x))[0]:x for x in session.odbs.keys()}
    max_loads = {}
    for job_name, odb_path in odb_names.items():
        print('Processing: {}'.format(job_name))
        odb = session.odbs[odb_path]
        get_crack_legnths(odb, job_name, method='path')
        process_outputs(odb, job_name, cutoff_freq)
        max_loads[job_name] = currentMax( session.xyDataObjects[job_name])[-1][1]
    
    # MBT
    repeated_portion_job_name = find_common_string(list(odb_names.keys()))
    added_mbt = False
    job_name_max_load = [k for k in max_loads.keys() if max_loads[k] == max(max_loads.values())]
    if len(job_name_max_load):
        job_name_max_load = job_name_max_load[0]
    else:
        job_name_max_load = ''
    for job_name, odb_path in odb_names.items():
        directory = os.path.dirname(odb_path)
        # Try to find matching par file
        par_file = [x for x in os.listdir(directory) if x == job_name + '.par']
        if len(par_file) == 1:
            par_file = os.path.join(directory, par_file[0])
            print('Using {} for MBT calculation'.format(par_file))
            param = get_parameters(par_file)
            param['af'] = param['a0']*2.2
            param['name'] = job_name.replace(repeated_portion_job_name, 'MBT_')
            if job_name == job_name_max_load:
                param['linear_segment'] = 'from_analysis'
                param['job_name'] = job_name
            mbt(**param)
            added_mbt = True

    # if not added_mbt:
    #     param = get_parameters('_DCB_3D_parameters.inp')
    #     param.update(get_parameters('_IM7_8552.inp'))
    #     param['af'] = param['a0']*2
    #     mbt(**param)
