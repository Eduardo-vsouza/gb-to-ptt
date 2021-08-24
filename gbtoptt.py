from Bio import SeqIO
import pandas as pd


class GenbankToPTT(object):
    def __init__(self, genbank_file):
        self.gb_record = SeqIO.read(open(genbank_file, "r"), "genbank")
        self.features = self.gb_record.features
        self.proteins = {'Location': [], 'Strand': [], 'Length': [], 'PID': [], 'Gene': [], 'Synonym': [], 'Code': [],
                     'COG': [], 'Product': []}
        self.rnas = {'Location': [], 'Strand': [], 'Length': [], 'PID': [], 'Gene': [], 'Synonym': [], 'Code': [],
                     'COG': [], 'Product': []}

    def get_records(self):
        for feature in self.features:
            # print('-----------------------------------------')
            # print(feature)
            # print(protein)
            # print(feature.qualifiers)
            # print(vars(feature))
            # print(feature.location)
            # print(feature.qualifiers.gene)
            # print(feature.strand)
            if feature.type == 'CDS':
                loc = f'{feature.location}'[1:-4].replace(":", "..")
                self.proteins['Location'].append(loc)
                self.proteins['Strand'].append(feature.strand)
                # print(vars(feature))
                g = False
                ltag = False
                prdct = False
                trslt = False
                for i in feature.qualifiers:
                    # if 'gene' not in feature.qualifiers:
                    #     print(feature.qualifiers)
                    # print(i)
                    feat = feature.qualifiers[i]
                    # print(feat)
                    if i == 'locus_tag':
                        ltag = True
                        self.proteins['Synonym'].append(feature.qualifiers[i][0])
                    if i == 'gene':
                        g = True
                        self.proteins['Gene'].append(feat[0])
                    if i =='product':
                        prdct = True
                        self.proteins['Product'].append(feat[0])
                    if i == 'translation':
                        trslt = True
                        self.proteins['Length'].append(len(feat[0]))
                        self.proteins["Code"].append("-")
                        self.proteins['PID'].append("-")
                        self.proteins['COG'].append("-")
                if g == False:
                    self.proteins['Gene'].append("-")
                if ltag == False:
                    self.proteins['Synonym'].append("-")
                if not prdct:
                    self.proteins['Product'].append("-")
                if not trslt:
                    self.proteins['Length'].append("-")
                    self.proteins["Code"].append("-")
                    self.proteins['PID'].append("-")
                    self.proteins['COG'].append("-")
                        # if i == 'note':
                        #     # print(feature.qualifiers[i])
                        #     loc = f'{feature.location}'[1:-4].replace(":", "..")
                        #     splat = loc.split("..")
                        #     start = int(splat[0])
                        #     end = int(splat[1])
                        #     length = end - start
                        #     # print(length)
                        #     self.proteins['Location'].append(loc)
                        #     # print(feature.location)
                        #     # print(feature.location)
                        #     self.proteins['Strand'].append(feature.strand)
                        #     # print(feature.qualifiers[i])
                        #     gene = feature.qualifiers[i][0].split(",")[0]
                        #     synonym = feature.qualifiers[i][0].split(",")[1]
                        #     self.proteins['Length'].append(length)
                        #     self.proteins['Gene'].append(gene)
                        #     self.proteins['Synonym'].append(synonym)
                        #     # self.rnas['Length'].append(feature.qualifiers[i][0].split(",")[3])
                        #     self.proteins["Code"].append("-")
                        #     self.proteins['PID'].append("-")
                        #     self.proteins['COG'].append("-")
                        #     self.proteins['Product'].append(gene)

                        # print(i)
                        # print(feature.qualifiers[i])
            else:
                for i in feature.qualifiers:
                    add = False
                    if i == 'note':
                        # print(feature.qualifiers[i])
                        if 'tRNA' in feature.qualifiers[i][0] or 'rRNA' in feature.qualifiers[i][0]:
                            loc = f'{feature.location}'[1:-4].replace(":", "..")
                            splat = loc.split("..")
                            start = int(splat[0])
                            end = int(splat[1])
                            length = end - start
                            # print(length)
                            self.rnas['Location'].append(loc)
                            # print(feature.location)
                            # print(feature.location)
                            self.rnas['Strand'].append(feature.strand)
                            # print(feature.qualifiers[i])
                            gene = feature.qualifiers[i][0].split(",")[0]
                            synonym = feature.qualifiers[i][0].split(",")[1]
                            self.rnas['Length'].append(length)
                            self.rnas['Gene'].append(gene)
                            self.rnas['Synonym'].append(synonym)
                            # self.rnas['Length'].append(feature.qualifiers[i][0].split(",")[3])
                            self.rnas["Code"].append("-")
                            self.rnas['PID'].append("-")
                            self.rnas['COG'].append("-")
                            self.rnas['Product'].append(gene)
                            # print(feature.qualifiers[note])

                    # if i == 'Product':
                    #     if add:
                    #         print(feature.qualifiers[i])
                    #         self.rnas["Product"].append(feature.qualifiers[i])
            # print(self.rnas)
            # print(self.proteins)

    def _add_features(self, feature, i, type='rna'):
        """

        :param feature: genbank feature
        :param type: either 'rna' or 'protein'
        :return:
        """
        ...
    def save_tables(self, out_ptt, out_rnt):

        ptt = pd.DataFrame(data=self.proteins)
        rnt = pd.DataFrame(data=self.rnas)
        ptt.to_csv(out_ptt, sep='\t', index=False)
        rnt.to_csv(out_rnt, sep='\t', index=False)


if __name__ == '__main__':
    data = GenbankToPTT(genbank_file='mtb.gb')
    data.get_records()
    data.save_tables(out_ptt='mtb.ptt', out_rnt='mtb.rnt')