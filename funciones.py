def get_eff_signal(self,datas,fs=512):
    stime = [x.get('Start_Time') for x in self.metadatas]
    etime = [x.get('End_Time') for x in self.metadatas]
    teff_s = np.max(stime)
    teff_e = np.min(etime)
    ts_i = [(teff_s - x).seconds for x in stime]
    tf_i = [(x - teff_e).seconds for x in etime]
    t_eff = (teff_e - teff_s).seconds
    data_eff = [x[fs * y:fs * (y + t_eff)] for x, y in zip(datas, ts_i)]
    data_eff = np.array(data_eff)
    print(t_eff)
    dt = 1 / fs
    tv = np.arange(0, t_eff, dt)
    return data_eff, t_eff, tv


def manage_data(self, paths):
    for path in paths:
        for path in paths:
            metadata = []
            with open(path, 'r', encoding='ISO-8859-1') as file:
                c = 0
                for line in file:
                    if line[0:6] != '[mm/s]':
                        tup = line.partition(':')

                        metadata.append([x.strip() for x in tup])

                        # print(line,c)
                        c += 1
                    else:

                        # print(line, c)
                        tup = line.partition(':')
                        metadata.append([x.strip() for x in tup])
                        c += 1
                        break
                # print('c', c)
    da, mda = create_data(path, metadata, c)
    datas.append(da)
    metadatas.append(mda)
    datas_eff, t_eff, tt = get_eff_signal()

def create_data(self, path, metadata, c, max=False):
        dic = {}
        dic.update({'Site_ID': metadata[0][2].split(',')[0]})
        dic.update({'Med_and_Node': metadata[0][2].split(',')[1].strip()[:-3]})
        dic.update({'Instrument': metadata[0][2].split(',')[1].strip()[-2:]})
        dic.update({'Instrument_Tag': metadata[2][2]})
        dic.update({'Node_Instrument': metadata[0][2].split(',')[1].split('_')[1] + '_' + metadata[0][2].split(',')[
                                                                                              1].strip()[-2:]})
        dic.update({'Data_Format_Bytes': int(metadata[4][2][:2])})
        dic.update({'Full_scale_mV': int(metadata[5][2][:3])})
        dic.update({'N_Channels': int(metadata[6][2])})
        dic.update({'Sampling_Rate': int(metadata[7][2][:3])})
        dic.update({'Start_Time': datetime.strptime(metadata[9][2], "%d/%m/%y %H:%M:%S")})
        dic.update({'End_Time': datetime.strptime(metadata[10][2], "%d/%m/%y %H:%M:%S")})
        dic.update({'Trace_length': metadata[11][2]})
        if dic['N_Channels'] == 4:
            print('channel 4')
            dic.update({'Start_Recording_UTC': metadata[18][2].split('\t')[1]})  ## will check this later
            dic.update({'Latitude': metadata[20][2]})
            dic.update({'Longitude': metadata[21][2]})
            dic.update({'Horizontal_Diluition': metadata[23][2]})
            dic.update({'Geoid_Altitude': metadata[24][2]})
            if not max:
                #print('q')
                data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2))
            if max:
                #print('qa')
                data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2), max_rows=max)
        if dic['N_Channels'] == 7:
            print('channel 7')
            dic.update({'Start_Recording_UTC': metadata[21][2].split('\t')[1]})  ## will check this later
            dic.update({'Latitude': metadata[23][2]})
            dic.update({'Longitude': metadata[24][2]})
            dic.update({'Horizontal_Diluition': metadata[26][2]})
            dic.update({'Geoid_Altitude': metadata[27][2]})
            if not max:
                data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2))
            if max:
                data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2), max_rows=max)
        if dic['N_Channels'] == 8:
            print('channel 8')
            dic.update({'Start_Recording_UTC': metadata[22][2].split('\t')[1]})  ## will check this later
            dic.update({'Latitude': metadata[24][2]})
            dic.update({'Longitude': metadata[25][2]})
            dic.update({'Horizontal_Diluition': metadata[27][2]})
            dic.update({'Geoid_Altitude': metadata[28][2]})
            if not max:
                data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2))
            if max:
                data = np.loadtxt(path, skiprows=c, encoding='ISO-8859-1', usecols=(2), max_rows=max)


        mdata_dic = dic
        return data, mdata_dic