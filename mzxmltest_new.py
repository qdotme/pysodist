from pyteomics import mzxml, pepxml
import numpy as np
import os
import time
import glob as glob
import matplotlib.pyplot as plt
MAX_SCANS = 2000  # maximum number of mzxml scans to look at once (must be greater than 50)
RETENTION_TIME_DELTA = 30  # retention time window to look in for summation integration
DELTA_MZ = 0.005  # adjusts the resolution of the mz array (filling in data points with interpolated values)
PEPTIDES = ["sp|P00359|G3P3_YEAST"]

def find_mzxml_pepxml():
    """
    function to find pepxml and mzxml files within the folder
    :return: tuple of two elements, first one which is the first mzxml file found and the second which is the first pep
    xml file found
    """
    mzxml_file_name = glob.glob("*.mz*")
    pepxml_file_name = glob.glob("*.pep.xml")
    return (mzxml_file_name,pepxml_file_name)


def list_summation_integration(scan_list, mass, charge, minrt, maxrt, peptide_name, dmz, label=None, min_mz=None, max_mz=None):
    """
    prints out a scatter plot of summation integration as per Davis's requirements
    :param scan_list: list of scans
    :param mass: mass found by the pepxml file
    :param charge: assumed charge from pepxml file
    :param minrt: minimum retention time to plot from
    :param maxrt: maximum retention time to plot to
    :param peptide_name: name of the peptide being graphed
    :param dmz: delta mz used to create y axis (default=.05)
    :param min_mz: minimum mz to graph from
    :param max_mz: maximum mz to graph to
    :param label: "SILAC-labeling"
    :return: prints out a scatter plot using parameters and Davis's requirements
    """

    if label == "silac-labeling":
#TODO: look at sequence, +8 lysine and +10 arginine
        silac_lower = int((mass+charge)/charge-5/charge) # always 5/charge for lower shift
        silac_upper = int((mass+charge)/charge+25/charge)
        return list_summation_integration(scan_list, mass, charge, minrt, maxrt, peptide_name, dmz, None, silac_lower, silac_upper)
    elif label == "15n-labeling":
#TODO: make a dictionary of amino acids and how many nitrogens it has, use to calculate m/z shift by looking at peptide sequence
        return list_summation_integration()
    elif label == "13c-labeling":
#TODO: same as above, but for carbon
        return list_summation_integration()

    mzarr = np.arange(float(min_mz),float(max_mz),dmz)
    mzlen = len(mzarr)
    intensity = np.zeros(mzlen)
    start_index = 0
    end_index = len(scan_list) - 1
    search_index = int(end_index/2)
    print("starting at rt:" + str(scan_list[search_index]['retentionTime']))
    while scan_list[search_index]['retentionTime'] > minrt:
        end_index = search_index
        search_index = int((start_index+end_index)/2)
    while scan_list[search_index]['retentionTime'] < minrt:
        start_index = search_index
        search_index = int((start_index+end_index)/2)
    while scan_list[search_index]['retentionTime'] > minrt:
        search_index -= 1

    i = int(search_index)

    print(scan_list[i]['retentionTime'],minrt)

    while scan_list[i]['retentionTime'] <= maxrt:
        print(scan_list[i]['retentionTime'])
        #interpIntensArr = np.interp(mzarr, scan_list[i]['m/z array'], scan_list[i]['intensity array'])
        #interp_intens_arr = interpolate_scan(scan_list[i], mass, charge, dmz, "silac-labeling")
        for x in range(0,mzlen-1):
            intensity[x] += scan_list[i]['intensity array'][x]
        i += 1
    print("ended at rt:" + str(scan_list[i]['retentionTime']))
    print(str(i-search_index)+" scans wide")

    filename = 'scans/test_data_' + str(peptide_name) + '_' + str(mass) + '.txt'

    plt.plot(mzarr,intensity)
    plt.show()
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename,'w+') as f:
        for j in range(len(intensity)):
            f.write(str(mzarr[j]) + " " + str(intensity[j]) + "\n")

    print("finished with:",str(peptide_name))


def interpolate_scan(scan, mass, charge, dmz, label=None, min_mz=None, max_mz=None):
    """
    returns a scan dictionary with its mz/intensity array interpolated to fit DELTA_MZ
    :param scan: dictionary object given by mzxml containing mz and intensity array
    :return: scan dictionary object with mz and intensity array interpolated to fit DELTA_MZ increments of mz
    """
    if label == "silac-labeling":
        silac_lower = int((mass+charge)/charge-5/charge)
        silac_upper = int((mass+charge)/charge+25/charge)
        return interpolate_scan(scan, mass, charge, dmz, None, silac_lower, silac_upper)
    elif label == "15n-labeling":
        return 1
    elif label == "13c-labeling":
        return 1
    mzarr = np.arange(float(min_mz), float(max_mz), dmz)
    interpIntensArr = np.interp(mzarr, scan['m/z array'], scan['intensity array'])
    scan['m/z array'] = mzarr
    scan['intensity array'] = interpIntensArr
    return scan

def find_first_hit(iter,minrt):
    """
    returns an iterator that starts with the first scan after the minimum retention_time
    :param iter: mzxmlfile iterator
    :param minrt: minimum retention_time (from first pepxml hit)
    :return: iterator at first scan after the minimum retention time
    """
    scan = iter.next()
    while scan['retentionTime']*60 < minrt:
        scan = iter.next()
    return iter


def run_program(mzxml_iterator,pepxml_iterator,scan_list=[]):
    try:
        pep_scan = pepxml_iterator.next()
    except:
        print("finished with pepxml")
        return time.time()

    #mzxml_iterator = find_first_hit(mzxml_iterator, pep_scan['retention_time_sec'] - RETENTION_TIME_DELTA)  # skip scans that will not be used.

    scan_list_len = len(scan_list)
    for i in range(MAX_SCANS-scan_list_len):
        try:
            scan = mzxml_iterator.next()
        except:
            print("finished mzxml file")
            return time.time()
        if i == MAX_SCANS - 50 - scan_list_len:
            ret_list = [scan]
        elif i > MAX_SCANS - 50 - scan_list_len:
            ret_list.append(scan)
        scan_list.append(scan)

    chosen_peptides = []  # [(peptide_sequence,retention_time,mz,z)]

    import csv
    with open('gmassacretrunc.input.csv', newline='\n') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        for row in spamreader:
            chosen_peptides.append((row[0],row[1],row[3],row[4]))
    j = 1
    while float(chosen_peptides[j][3]) * 60 < scan_list[-1]['retentionTime'] * 60 - 50:
        peptide_sequence = chosen_peptides[j][0]
        print("Peptide sequence: " + peptide_sequence)
        precursor_mass = float(chosen_peptides[j][2])*float(chosen_peptides[j][3])-float(chosen_peptides[j][3])
        assumed_charge = int(chosen_peptides[j][3])

        mz = float(chosen_peptides[j][2])
        print('mz from csv', mz)
        mz = round(mz, 3)
        t_mz = mz * 1000
        t_mz += 5 - t_mz % 5
        mz = t_mz / 1000  # is the mz that fits as a multiple of mz
        int_rt_arr = []
        rt_arr = []
        print('mz=',mz)

        # interpolate the scan_list
        for index, value in enumerate(scan_list):
            scan_list[index] = interpolate_scan(scan_list[index], precursor_mass, assumed_charge, DELTA_MZ,
                                                "silac-labeling")  # interpolate the scan
            index_found = int((mz - scan_list[index]['m/z array'][0]) / (DELTA_MZ))  # find the index of mz in the scan

            print(index_found)
            int_rt_arr.append(
                scan_list[index]['intensity array'][index_found])  # add the intensity at that mz to the array
            rt_arr.append(scan_list[index]['retentionTime'])  # add the rt to the list of rts

        plt.scatter(rt_arr, int_rt_arr)
        plt.show()
        # find critical points (peaks)
        peak_index = 0
        peak_value = 0
        peak_rt = 0
        peak_half_max_rt = 0
        peak_width = 0

        for index, intensity in enumerate(int_rt_arr):
            if intensity < int_rt_arr[index - 1] and (
                    int_rt_arr[index - 1] > peak_value):  # the peak went fron increasing to decreasing
                peak_index = index - 1
                peak_value = int_rt_arr[index - 1]
                peak_rt = rt_arr[index - 1]
            if peak_index > 0 and (int_rt_arr[index - 1] > peak_value / 2 > intensity):
                peak_half_rt = (rt_arr[index - 1] + rt_arr[index]) / 2
                # 1/4 width = peak_half_rt - peak_rt
                peak_width = 4 * (peak_half_rt - peak_rt)

        print("rt of peak: " + str(peak_rt), "full width of peak in rt:" + str(peak_width))

        start_rt = peak_rt - peak_width/2
        end_rt = peak_rt + peak_width/2

        print("start rt %s, end rt %s" % (start_rt, end_rt))

        list_summation_integration(scan_list, precursor_mass, int(chosen_peptides[j][3]),
                                   start_rt, end_rt, peptide_sequence, DELTA_MZ,
                                   "silac-labeling")  # HARD CODED SILAC LABELING

        data_file_name = 'scans/test_data_' + str(peptide_sequence) + '_' + str(precursor_mass) + '.txt'

        with open('misc/NPULSE.batch', 'a') as f:
            f.write(str(peptide_sequence) + " " + str(precursor_mass) + " " + str(
                data_file_name) + "\n")

        j += 1

#pep_scan['retention_time_sec'] < scan_list[-1]['retentionTime'] * 60 - 50
    while False:
        peptide_name = pep_scan['search_hit'][0]['proteins'][0]['protein']
        if "RAND" not in peptide_name and peptide_name in PEPTIDES:  # look for peptides

            print("Peptide: " + pep_scan['search_hit'][0]['peptide'])

            id_time = pep_scan['retention_time_sec']
            precursor_mass = pep_scan['precursor_neutral_mass']
            assumed_charge = pep_scan['assumed_charge']

            mz = (precursor_mass + assumed_charge) / assumed_charge
            mz = round(mz, 3) 
            t_mz = mz * 1000
            t_mz += 5 - t_mz % 5
            mz = t_mz / 1000  # is the mz that fits as a multiple of mz
            int_rt_arr = []
            rt_arr = []

            # interpolate the scan_list
            for index, value in enumerate(scan_list):
                scan_list[index] = interpolate_scan(scan_list[index], precursor_mass, assumed_charge, DELTA_MZ, "silac-labeling")  # interpolate the scan
                index_found = int((mz - scan_list[index]['m/z array'][0])/(DELTA_MZ))  # find the index of mz in the scan
                int_rt_arr.append(scan_list[index]['intensity array'][index_found])  # add the intensity at that mz to the array
                rt_arr.append(scan_list[index]['retentionTime'])  # add the rt to the list of rts

            plt.scatter(rt_arr, int_rt_arr)
            plt.show()
            # find critical points (peaks)
            peak_index = 0
            peak_value = 0
            peak_rt = 0
            peak_half_max_rt = 0
            peak_width = 0

            for index, intensity in enumerate(int_rt_arr):
                if intensity < int_rt_arr[index-1] and (int_rt_arr[index - 1] > peak_value): # the peak went fron increasing to decreasing
                    peak_index = index - 1
                    peak_value = int_rt_arr[index-1]
                    peak_rt = rt_arr[index-1]
                if peak_index > 0 and (int_rt_arr[index-1] > peak_value/2 > intensity):
                    peak_half_rt = (rt_arr[index-1] + rt_arr[index])/2
                    # 1/4 width = peak_half_rt - peak_rt
                    peak_width = 4 * (peak_half_rt - peak_rt)

            print(peak_rt,peak_width)

            start_rt = peak_rt - 3 * peak_width
            end_rt = peak_rt + 10 * peak_width

            print(start_rt,end_rt)

            mz_found = int(
                (pep_scan['precursor_neutral_mass'] + pep_scan['assumed_charge']) / pep_scan['assumed_charge'])

            list_summation_integration(scan_list, pep_scan['precursor_neutral_mass'], pep_scan['assumed_charge'],
                                       start_rt, end_rt, peptide_name, DELTA_MZ, "silac-labeling")  # HARD CODED SILAC LABELING

            data_file_name = 'scans/test_data_' + str(peptide_name) + '_' + str(pep_scan['precursor_neutral_mass']) + '.txt'

            with open('misc/NPULSE.batch', 'a') as f:
                f.write(str(pep_scan['search_hit'][0]['peptide']) + " " + str(pep_scan['assumed_charge']) + " " + str(data_file_name) + "\n")

        try:
            pep_scan = pepxml_iterator.next()
        except:
            print("finished with pepxml file")
            return time.time()
    del scan_list[:]
    del scan_list
    return run_program(mzxml_iterator,pepxml_iterator,ret_list)

results = find_mzxml_pepxml()
print("Found mzxml file:", results[0][0], "and pepxml file", results[1][0])

mzxmlfilename = results[0][0]
print("accessing mzxml file:", mzxmlfilename)
mzxml_it = mzxml.MzXML(mzxmlfilename)

pepxmlfilename = results[1][0]
print('accessing pepxml file:', pepxmlfilename)
try:
    pepxml_it = pepxml.PepXML(pepxmlfilename)
except:
    pepxml_it = pepxml.PepXML(pepxmlfilename)
os.makedirs(os.path.dirname('misc/NPULSE.batch'), exist_ok=True)  # create NPULSE file
t1 = time.time()
t2 = run_program(mzxml_it, pepxml_it)
print(t1,t2)
print(t2 - t1)
dt = t2 - t1
print("time taken:",str(dt))
