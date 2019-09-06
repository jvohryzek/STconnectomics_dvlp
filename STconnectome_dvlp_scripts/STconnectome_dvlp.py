#!/usr/bin/python

import os
import networkx as nx
import numpy as np
import scipy.io as mat
import sklearn as sk
from sklearn import preprocessing
import nibabel as nib
import csv
# ======================================================================

''' Ia. INPUT and OUTPUT directories'''

SC_filename = "STconnectome_dvlp_scripts/INPUT/SCgroup50_cortical.txt"
LABEL_filename = "STconnectome_dvlp_scripts/INPUT/brain_labels_448.txt"
cortical_indexes_file = "STconnectome_dvlp_scripts/INPUT/index_CORTICAL_Laus2008_scale4.mat"
list_filename = "STconnectome_dvlp_scripts/INPUT/subjects_list.txt"
yeo_filename = "STconnectome_dvlp_scripts/INPUT/Lausanne2008_Yeo7RSN_scale4"
age_filename = "STconnectome_dvlp_scripts/INPUT/Age_subject_age.txt"
age_name_filename = "STconnectome_dvlp_scripts/INPUT/Age_subject_names.txt"

''' Ib. INPUT and OUTPUT directories'''
output_dir_feature_matrix = "STconnectome_dvlp_scripts/OUTPUT/feature_matrix"
output_dir_feature_matrix_filtered = "STconnectome_dvlp_scripts/OUTPUT/feature_matrix_filtered"
output_dir_cc_filtered = "STconnectome_dvlp_scripts/OUTPUT/cc_filtered"

''' II. Indices of cortical regions'''
# excluding subcortical regions
ix_cort = mat.loadmat(cortical_indexes_file)
ix_cort = np.squeeze(np.add(ix_cort.get('ix'), -1))

'''III. Yeo atlas'''
yeo = open(yeo_filename, 'r').read().split('\n')
yeo = [yeo[i] for i in ix_cort]

''' IVa. SC Data !!!WORKING WITH!!!'''
SC = np.loadtxt(open(SC_filename, "rb"), dtype=int, delimiter=',', skiprows=0)
# Remove diagonal elements (not meaningful for the construction of a spatio-temporal connectome)
SC = SC - np.diag(np.diag(SC))
''' V. Subjects age'''
names = []
with open(age_name_filename, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        names.append(row[0])
age = []
with open(age_filename, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        age.append(row[0])

''' VI. Initiate Parameters'''
ts = []  # timeseries
patients = []
count = 0
pos_threshold = 2.03

'''START OF THE LOOP'''
# set path
base_input_path = '/Volumes/JakubExtHD/RAD/ELLIEZ_DEVELOPMENT/DATA'
# browse input folder
files_list = os.listdir(base_input_path)
count = 0
count_count = 0

for root, dirs, files in os.walk(base_input_path):
    for filename in dirs:
        if filename.find('3T') != -1 & filename.find('5868_3_3T') == -1 & filename.find('5227_3_3T') == -1:  # select only file with name 'sbj*' in input folder / excluding elliez and the wierd NAN in time series subject for now

            count = count+1 # counts the subjects for the original data set of 101

            print('Subject ' + str(count) + ': ' + filename + ' , Age: ' + age[names.index(filename[0:6])])
            if os.path.exists(os.path.join(base_input_path, filename, 'T1/CMP/fMRI/preprocessing_Jakub_v3/averageTimeseries_250_FIR_filtered.mat')): # before the excluded sbj but when do home scrubbing we take only the sbj after exlusion so CCfilt)scrub has the right amount of sbj
                path = os.path.join(base_input_path, filename, 'T1/CMP/fMRI/preprocessing_Jakub_v3/averageTimeseries_250_FIR_filtered.mat')
                path_ICV = os.path.join(base_input_path, filename, 'T1/FREESURFER/stats')
                count_count = count_count + 1 # counts the subjects after scrubbing
                print(count_count)
                if os.path.isfile(path):
                    ''' 1. load timeseries data '''
                    ts_load = mat.loadmat(path)
                    ts = ts_load['ts_fir']
                    ts = ts[np.ix_(ix_cort, np.arange(0, ts.shape[1]))]  # just gets rid of the subcortical regions resulting in ndarray of 448x196

                    subject = filename[0:10]
                    patients.append(filename[0:10])

                    ''' 2. ROI labels '''
                    labels = open(LABEL_filename, 'r').read().split('\n')
                    if not labels[-1]:
                        labels = labels[0:-1]

                    # Number of ROIs in the structural connectivity graph
                    nROIs = ts.shape[0]
                    # Number of time points
                    ntp = ts.shape[1]

                    '''3. Intracranial Volume (ICV)'''

                    # I. version

                    # Extract ICV information
                    data = []
                    with open(os.path.join(path_ICV, 'aseg.stats'), 'r') as f:
                        for line in f:
                            if line.__contains__('IntraCranialVol'):
                                start = line.find('Volume, ') + 8
                                end = line.find(', mm')
                                ICV = float(line[start:end])

                    # III. version

                    data_aseg = nib.load(os.path.join(base_input_path, filename, 'T1/FREESURFER/mri/aseg.nii.gz'))
                    mask_aseg = data_aseg.get_data()
                    ICV_aseg = sum(sum(sum(mask_aseg > 0))) # in mm^3


                    ''' 4. Define output file names'''
                    # ======================================================================

                    FM_filename = os.path.join(output_dir_feature_matrix, subject + "_FM.mat")
                    FMfilt_filename = os.path.join(output_dir_feature_matrix_filtered, subject + "_FMfilt.mat")
                    CCfilt_filename = os.path.join(output_dir_cc_filtered, subject + "_CCfilt.mat")

                    ''' 5. z-score ROI-wise time series '''
                    # ======================================================================
                    # print('..... z-score time series .....')
                    test = np.where(ts == 0)
                    test = np.array(test[0])
                    testn = test.size
                    test = np.where(ts == 1)
                    test = np.array(test[0])
                    testn = testn + test.size
                    if ts.size != testn:
                        for i in range(len(ts)):
                            ts[i] = sk.preprocessing.scale(ts[i], with_mean = True, with_std = True, copy = True)
                    else:
                        pos_threshold = 1

                    # print('..... write file: ' + TS_zscore_filename)

                    ''' 6. Create static structural connectivity graph Gs from the input adjacency matrix SC '''
                    # ======================================================================
                    # Create the base graph from the adjacency matrix
                    Gs = nx.from_numpy_matrix(SC)

                    ''' 7. Build a spatio-temporal connectome from SC_nof and ts information '''
                    # ======================================================================
                    # print('..... build multilayer network (spatio-temporal connectome) .....')
                    # Initiate a new NetworkX graph
                    G = nx.Graph()
                    graph_data = {}
                    graph_data['subject'] = subject  # subject IDs
                    graph_data['ntp'] = ntp  # ntp, number of time points
                    graph_data['nROIs'] = nROIs  # nROIs, number of anatomical ROIs
                    graph_data['activation_threshold'] = pos_threshold  # nROIs, number of anatomical ROIs
                    graph_data['patient_age'] = age[names.index(filename[0:6])]  # age of the subjects
                    graph_data['Intracranial_Volume'] = ICV # in mm^3
                    graph_data['Intracranial_Volume_aseg'] = ICV_aseg  # in mm^3
                    G.graph['graph_data'] = graph_data

                    # Loop throughout all the time points (1 -> ntp-1)
                    for t in range(ntp):

                        # For current time point, find active ROIs
                        tsst = ts[:, t]  # fMRI values for all the nodes, current time point (t)
                        active_nodes = np.where(tsst >= pos_threshold)[0]

                        # Loop throughout all the ROIs active at current time point
                        for i in active_nodes:

                            # Generate a new node ID (in the multilayer network)
                            # NOTE: each node in the multilayer network has a unique ID equal to layer_pos * nROIs + i,
                            # with i anatomical_id of the considered node (from 1 to nROIs), layer_pos node position in
                            # time (from 1 to ntp), and nROIs number of ROIs in the structural connectivity graph
                            node_id = t * nROIs + i + 1  # ROI IDs start from 1
                            # If node_id does not exist in G, add it
                            if ~G.has_node(node_id):
                                # Generate attributes for the new node in the multilayer network
                                node_attrs = {}
                                node_attrs['anatomical_id'] = i + 1  # ROIs ID start from 1
                                node_attrs['weight'] = tsst[i]
                                node_attrs['tp'] = t + 1  # tp ID start from 1
                                node_attrs['node_id'] = node_id
                                node_attrs['label'] = labels[i]
                                node_attrs['patient_age'] = age[names.index(filename[0:6])]
                                node_attrs['functional_network'] = yeo[i]  # functional network according to yeo 7 f. net
                                node_attrs['Intracranial_Volume'] = ICV  # in mm^3
                                node_attrs['Intracranial_Volume_aseg'] = ICV_aseg  # in mm^3
                                G.add_node(node_id, node_attrs)

                            if t < (ntp - 1):
                                # Extract the neighbors of anatomical_id in Gs (structural connectivity graph)
                                neighbor_nodes = Gs.neighbors(i)
                                # Consider as well the node itself in the following (t+1) time point
                                neighbor_nodes.extend([i])

                                # Extract brain regions that are active at the following (t+1) time point
                                tsstt = ts[:, t + 1]  # fMRI values for all the nodes, following (t+1) time point
                                active_nodes_tt = np.where(tsstt >= pos_threshold)[0]

                                '''Important line'''
                                # Intersect current region's neighbors, and regions active at the following (t+1) time point
                                new_nodes = np.intersect1d(neighbor_nodes, active_nodes_tt, assume_unique=True)

                                # Loop throughout all neighbor and active regions: add them to the multilayer network, together with the corresponding edges
                                for j in new_nodes:
                                    node_id_new = ((t + 1) * nROIs) + j + 1
                                    node_attrs = {}
                                    node_attrs['anatomical_id'] = j + 1
                                    node_attrs['weight'] = tsstt[j]
                                    node_attrs['tp'] = t + 2
                                    node_attrs['node_id'] = node_id_new
                                    node_attrs['label'] = labels[j]
                                    node_attrs['patient_age'] = age[names.index(filename[0:6])]
                                    node_attrs['functional_network'] = yeo[j] # functional network according to yeo 7 f. net
                                    node_attrs['Intracranial_Volume'] = ICV  # in c^3
                                    node_attrs['Intracranial_Volume_aseg'] = ICV_aseg  # in mm^3
                                    G.add_node(node_id_new, node_attrs)
                                    G.add_edge(node_id, node_id_new)

                    ''' 8. Add bi-directional links between nodes belonging to the same layer of Gs'''
                    # TO NOTE HERE: the bid-directional links are encoded so that we have two edges say (node1, node2) and (node2,node1)
                    # Loop throughout all the time points (1 -> ntp)
                    for t in range(ntp):

                        # For current time point, find active nodes
                        tsst = ts[:, t]  # fMRI values for all the nodes, current time point
                        active_nodes = np.where(tsst >= pos_threshold)[0]

                        # Add links between active nodes at current time point, which are also neighbors in Gs
                        for i in active_nodes:

                            # Node ID in multilayer network
                            node_id = t * nROIs + i + 1

                            # Extract region neighbors in Gs
                            neighbor_nodes = Gs.neighbors(i)
                            '''Important line'''
                            # Intersect current node neighbors, and nodes active at the current (t) time point
                            new_nodes = np.intersect1d(neighbor_nodes, active_nodes, assume_unique=True)

                            # Add edges
                            for j in new_nodes:
                                node_id_new = t * nROIs + j + 1
                                if ~G.has_edge(node_id, node_id_new):
                                    G.add_edge(node_id, node_id_new)

                    ''' 9. Save spatio-temporal connectome (multilayer network) as gpickle file '''
                    # print('..... write file: ' + G_filename)
                    # nx.write_gpickle(G, G_filename)

                    ''' 10. Extract the connected components (CCs) of the multilayer network '''
                    # ======================================================================
                    # print('..... extract connected components (CCs) of multilayer network .....')
                    # The output nx.connected_component_subgraphs() is a list of nx graphs,
                    # ordered from the largest to the smallest (in terms of number of nodes)
                    CC = list(nx.connected_component_subgraphs(G))
                    # print('.....    (number of CCs: ' + str(len(CC)) + ')')

                    # Set attributes of CCs: width (temporal extension), height (spatial extension), subject ID
                    # And save the CCs to a Matlab array of structures
                    height = np.zeros(shape=(len(CC)))
                    width = np.zeros(shape=(len(CC)))
                    subjectCC = []
                    subjectCC_age = []

                    # Initialize Python list of dictionaries and feature matrix
                    dictlist_cc = [dict() for x in range(len(CC))]
                    FM = np.zeros(shape=(len(CC), nROIs))

                    # Loop throughout all the connected components
                    for i in range(0, len(CC)):

                        # Current connected component
                        cc = CC[i]
                        # Nodes anatomical_id and layer_pos
                        nodes = cc.nodes()

                        functional_network_dict = nx.get_node_attributes(cc, 'functional_network')
                        functional_network = [functional_network_dict[x] for x in nodes]

                        anatomical_id_dict = nx.get_node_attributes(cc, 'anatomical_id')
                        anatomical_id = [anatomical_id_dict[x] for x in nodes]

                        tp_dict = nx.get_node_attributes(cc, 'tp')
                        tp = [tp_dict[x] for x in nodes]

                        # Spatial height and temporal width of current component
                        height[i] = len(set(anatomical_id))
                        width[i] = max(tp) - min(tp) + 1
                        # Subject
                        # subjectCC.append(subject)
                        # Feature vector
                        ids = np.array(anatomical_id) - 1
                        for j in range(0, len(ids)):
                            FM[i][ids[j]] += 1

                        # Fill-in Python dictionary
                        dictlist_cc[i] = {'height': height[i], 'width': width[i], 'subject': subject,
                                          'subject_age': age[names.index(filename[0:6])], 'nodes': nodes, 'edges': cc.edges(),
                                          'anatomical_id': anatomical_id,'functional_network': functional_network,
                                          'tp': tp, 'ICV':ICV,'ICV_aseg':ICV_aseg}

                    # Save connected components list of dictionaries to Matlab format
                    # print('..... write file: ' + CC_filename)
                    mdict = {'CC': dictlist_cc}
                    # mat.savemat(CC_filename, mdict)

                    # Normalize and save feature matrix
                    # print('..... write file: ' + FM_filename)
                    FM_norms = np.apply_along_axis(np.linalg.norm, 1, FM)
                    FM = FM.astype(float) / FM_norms[:, np.newaxis]
                    mat.savemat(FM_filename, mdict={'FM': FM})

                    # Filter connected components according to their height and width
                    # [THIS PART WILL BE REMOVED IN THE PUBLIC DISTRIBUTION OF THE SCRIPT!]
                    # ======================================================================

                    filt_min_width = 2  # WIDTH is the temporal span of the CC (number of time points)
                    filt_max_width = 276
                    filt_min_height = 6  # HEIGHT is the spatial span of the CC (number of anatomical regions)
                    filt_max_height = 448
                    filt_min_nn = 6  # overall number of nodes in CC
                    filt_max_nn = 276 * 448

                    CCfilt = []
                    for i in range(0, len(CC)):
                        cc = CC[i]
                        if filt_min_width <= width[i] <= filt_max_width:
                            if filt_min_height <= height[i] <= filt_max_height:
                                if filt_min_nn <= cc.number_of_nodes() <= filt_max_nn:
                                    CCfilt.append(CC[i])
                    # print('..... (number of filtered CCs: ' + str(len(CCfilt)) + ')')

                    # Initialize Python list of dictionary
                    dictlist_ccfilt = [dict() for x in range(len(CCfilt))]
                    FMfilt = np.zeros(shape=(len(CCfilt), nROIs))
                    height_filt = np.zeros(shape=(len(CCfilt)))
                    width_filt = np.zeros(shape=(len(CCfilt)))
                    subjectCC_filt = []
                    subjectCC_filt_age = []
                    # Loop throughout all the connected components
                    for i in range(0, len(CCfilt)):

                        # Current connected component
                        cc = CCfilt[i]
                        # Nodes anatomical_id and layer_pos
                        nodes = cc.nodes()

                        functional_network_dict = nx.get_node_attributes(cc, 'functional_network')
                        functional_network= [functional_network_dict[x] for x in nodes]

                        anatomical_id_dict = nx.get_node_attributes(cc, 'anatomical_id')
                        anatomical_id = [anatomical_id_dict[x] for x in nodes]

                        tp_dict = nx.get_node_attributes(cc, 'tp')
                        tp = [tp_dict[x] for x in nodes]

                        # Spatial height and temporal width of current component
                        height_filt[i] = len(set(anatomical_id))
                        width_filt[i] = max(tp) - min(tp) + 1

                        # Feature vector
                        ids = np.array(anatomical_id) - 1
                        for j in range(0, len(ids)):
                            FMfilt[i][ids[j]] += 1

                        # Fill-in Python dictionary
                        dictlist_ccfilt[i] = {'height': height_filt[i], 'width': width_filt[i],
                                              'subject': subject, 'subject_age': age[names.index(filename[0:6])],
                                              'nodes': nodes, 'edges': cc.edges(), 'anatomical_id': anatomical_id,
                                              'functional_network': functional_network,'tp': tp, 'ICV': ICV, 'ICV_aseg': ICV_aseg}

                    # Save connected components list of dictionaries to Matlab format
                    print('..... write file: ' + CCfilt_filename)
                    mdict = {'CCfilt': dictlist_ccfilt}
                    mat.savemat(CCfilt_filename, mdict)
                    # Normalize and save feature matrix
                    print('..... write file: ' + FMfilt_filename)
                    FMfilt_norms = np.apply_along_axis(np.linalg.norm, 1, FMfilt)
                    FMfilt = FMfilt.astype(float) / FMfilt_norms[:, np.newaxis]
                    mat.savemat(FMfilt_filename, mdict={'FMfilt': FMfilt})

                    # Clean
                    del G
                    del ts
                    del height
                    del height_filt
                    del width
                    del width_filt
