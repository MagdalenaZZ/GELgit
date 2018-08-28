

/genomes/software/src/test-venvs/cancer-test/local/lib/python2.7/site-packages/CancerPlotGeneration/vcf_tools.py


/genomes/software/src/test-venvs/cancer-test/local/lib/python2.7/site-packages/DataModels/gVCFVariant.py

import copy


__author__ = 'antonior'


class GelVariant(Variant):
    def __init__(self, filter=None, qual=0, annotation=None, **kwargs):
        """

        :type annotation: GelVariantAnnotation
        """
        self.annotation = annotation
        self.qual = qual
        self.filter = list(filter) if filter is not None else ['PASS']
        super(GelVariant, self).__init__(**kwargs)
        self.type = self.define_type
        self.sample_names = [call.callSetName for call in self.calls]

        

    def get_sample_value(self, sample_name, field_name):
        """

        Get the specified field value for the sample given in that variant

        :type field_name: str
        :type sample_name: str
        :rtype : str
        :param sample_name: Query sample
        :param field_name: The name of Field in Format will be returned if exists (i,e DP4)
        :return: Information stored in a format field given in a sample given
        """
        sample = self.get_call(sample_name)
        if field_name in sample.info:
            return sample.info[field_name]
        else:
            raise StandardError('The sample field "' + field_name +
                                '" is not present in this variant')



    def calculate_vaf(self, sample_name):
        vaf = dict()
        if self.type != "Indel":
            try:
                if isinstance(self.alternateBases, list):
                    alt = self.alternateBases[0]
                else:
                    alt = self.alternateBases

                ref = self.referenceBases

                try:
                    AU = self.get_sample_value(sample_name, "AU")
                except:
                    AU = ['0', '0']
                try:
                    CU = self.get_sample_value(sample_name, "CU")
                except:
                    CU = ['0', '0']
                try:
                    GU = self.get_sample_value(sample_name, "GU")
                except:
                    GU = ['0', '0']
                try:
                    TU = self.get_sample_value(sample_name, "TU")
                except:
                    TU = ['0', '0']
                # this uses all reads at that position (more risky but more accurate?)
                if "." not in AU:
                    a = int(AU[0])
                else:
                    a = 0
                if "." not in CU:
                    c = int(CU[0])
                else:
                    c = 0
                if "." not in GU:
                    g = int(GU[0])
                else:
                    g = 0
                if "." not in TU:
                    t = int(TU[0])
                else:
                    t = 0
                counts = {
                    "A": a,
                    "C": c,
                    "G": g,
                    "T": t
                }
                # DP = int(self.get_sample_value(sample_name, "DP")[0])
                # try:
                #     FDP = int(self.get_sample_value(sample_name, "FDP")[0])
                # except:
                #     print self

                var_count = int(counts[alt])
                total_count = int(counts[alt]) + int(counts[ref])
                vaf["ref"] = total_count - var_count
                vaf["alt"] = var_count
                if total_count > 0:
                    vaf["fraction"] = (var_count / float(total_count))
                else:
                    vaf["fraction"] = 0
            except:
                vaf["ref"] = -1
                vaf["alt"] = -1
                vaf["fraction"] = -1
        else:
            TIR = int(self.get_sample_value(sample_name, "TIR")[0])
            TAR = int(self.get_sample_value(sample_name, "TAR")[0])
            vaf["ref"] = TAR
            vaf["alt"] = TIR
            try:
                vaf["fraction"] = (TIR / float(TAR + TIR))
            except:
                vaf["fraction"] = 0
        return vaf





