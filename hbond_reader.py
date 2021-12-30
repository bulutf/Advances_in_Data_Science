import pandas as pd
import re

class hbond_reader:
        
    def write_avghbond(self):

        df = pd.read_table(self.filename, sep="\s+")
        transpose=df.T
        transpose=transpose.drop('#Frame')
        
        seperated = re.findall("(\w+_\d+@\w+)-(\w+_\d+@\w+)-(\w+)",str(list(transpose.index)))

        acceptor=list()
        donor=list()
        donorh=list()

        for i in range(len(transpose.index)):
            acceptor.append(seperated[i][0])
            donor.append(seperated[i][1])
            donorh.append(seperated[i][1].replace('@N', '@HN'))

        frames=list()
        frac=list()

        for i in range(len(transpose.sum(axis=1))):
            frames.append(transpose.sum(axis=1)[i])
            frac.append(transpose.sum(axis=1)[i]/200)

        d = {'#Acceptor': acceptor, 'DonorH': donorh, 'Donor':donor, 'Frames':frames,'Frac':frac}
        df = pd.DataFrame(data=d)
        df.to_csv('avghbond.dat', index=False)

    def write_hbondnumber(self):
        
        df = pd.read_table(self.filename, sep="\s+")
        transpose=df.T
        transpose=transpose.drop('#Frame')
        
        frames=list(df["#Frame"])
        frames
            
        tseries=list()
        for i in range(len(transpose.columns)):
            tseries.append(transpose.sum(axis=0)[i])
            
        d = {'Frames': frames, 'Total_Number_of_Hbonds': tseries}
        df = pd.DataFrame(data=d)
        df.to_csv('hbond_number_tseries.dat', index=False)
        
    def check_intraaminoacid(self):
        
        df = pd.read_table(self.filename, sep="\s+")
        transpose=df.T
        transpose=transpose.drop('#Frame')
        seperated2 = re.findall("\w+_(\d+)@\w+-\w+_(\d+)@\w+-\w+",str(list(transpose.index)))
        
        hbonds=list()
        numberofframes=list()

        for i in range(len(seperated2)):
            if seperated2[i][0]==seperated2[i][1]:
                hbonds.append(transpose.index[i])
                numberofframes.append(transpose.sum(axis=1)[i])

        if len(hbonds)==0:
            print("no intra amino acid hbond found")
        else:
            d = {'hbonds': hbonds, '# of frames': numberofframes}
            df = pd.DataFrame(data=d)
            print(df)
            
    def find_highest_donortype(self):
        
        df = pd.read_table(self.filename, sep="\s+")
        transpose=df.T
        transpose=transpose.drop('#Frame')
        
        seperated = re.findall("(\w+_\d+@\w+)-(\w+_\d+@\w+)-(\w+)",str(list(transpose.index)))

        donor=list()

        for i in range(len(transpose.index)):
            donor.append(seperated[i][1])

        print("Amino acid type that gives the highest number of hbond donor:\n{}".format(max(set(donor), key = donor.count)))
