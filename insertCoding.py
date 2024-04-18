# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 11:16:45 2024

@author: Wei Qiang
"""
import sys, copy, math, array, io, os, cv2, random
from itertools import combinations
from md5hash import scan
from PIL import Image
from pydub import AudioSegment
import numpy as np
import scipy.io.wavfile
import moviepy.editor as mp
import imageio
from multiprocessing import Pool
from multiprocessing import Lock
import itertools
import argparse


class tool:

    def decimal2OtherSystem(self, decimal_number: int, other_system:int, precision:int = None) -> list:
        '''
        list中，每一个元素为一位，例十六进制中，"10"="A"
        precision: 转换的位数，也就是序列长度，如果有指定，则在前面补0
        '''

        remainder_list = []
        while True:
            remainder = decimal_number%other_system
            quotient = decimal_number//other_system
            remainder_list.append(str(remainder))
            if quotient == 0:
                break
            decimal_number = quotient

        # if other_system <10:
        #     other_system_num = "".join(remainder_list[::-1])
        #     return other_system_num
        # else:
            # return remainder_list[::-1]

        num_list = remainder_list[::-1]
        # 指定precision
        if precision != None:
            if precision < len(num_list):
                raise ValueError("The precision is smaller than the length of number. Please check the [precision] value!")
            else:
                num_list = ["0"]*(precision - len(num_list)) + num_list
        return num_list


    def convert2Decimal(self, num, system: int) -> int:
        num_len = len(num)
        dec_num = 0
        for i in range(num_len):
            dec_num += int(num[i])*system**(num_len -(i+1))

        return dec_num


    def getLCM(self, input_system: int, convert_system: int):

        file_bits = math.log(input_system, 2)
        ## 当需要转换的system不是2的指数时，取整所得的最小公倍数，可让实际转换后单位数量中可包含的信息量大于文件的信息量
        sys_bits = int(math.log(convert_system, 2))

        n = 1
        while True:
            lowest_common_multiple = file_bits*n
            if lowest_common_multiple%sys_bits ==0:
                break
            else:
                n += 1

        file_bits_num = int(lowest_common_multiple/file_bits)
        sys_bits_num = int(lowest_common_multiple/sys_bits)

        return file_bits_num, sys_bits_num


    def xor(self, seq_list, digit = 4):
        len_list = set([len(i) for i in seq_list])
        if len(len_list) != 1:
            raise ValueError("The length of two sequences is not same!!\n")

        xor_seq = []
        for i in range(len(seq_list[0])):
            xor_num = sum([int(j[i]) for j in seq_list])%digit
            xor_seq.append(str(xor_num))


        return xor_seq

    def xorDecode(self, xor_seq, seq_1, digit = 4):
        if len(xor_seq) != seq_1:
            raise ValueError("The length of two sequences is not same!!\n")

        seq_2_list = []
        for i in range(xor_seq):
            if int(xor_seq[i]) > int(seq_1[i]):
                seq_2_list.append(str(int(xor_seq[i]) - int(seq_1[i])))
            elif int(xor_seq[i]) < int(seq_1[i]):
                seq_2_list.append(str(digit - seq_1[i] + xor_seq[i]))

        return seq_2_list


    def getTxtCodingDic(self, system):
        # 目前生成4组映射数值，一般是够用的,每一组根据system选择字符个数
        total_num = 0
        bytes_num_dic = {}

        if system >= 256:
            cha = 1
            for i in range(4):
                max_num = [system-2]*(i+1)
                max_num = self.convert2Decimal(max_num, system)
                total_num += max_num+1
                bytes_num_dic[total_num] = i+1

        else:
            if system == 4:
                cha = 4
            elif system == 8:
                cha = 3
            elif system in [16, 32, 64, 128]:
                cha = 2
            
            for i in range(4):
                # max_num is list
                max_num = [system-2] + (cha-1)*(i+1)*[system - 1]
                max_num = self.convert2Decimal(max_num, system)
                total_num += max_num+1
                bytes_num_dic[total_num] = i+1


        return bytes_num_dic, cha

class insertEncode:

    def __init__(self, input_path, output_dir, split_len=200, homo=4, gc_standard = 0.5, redanduncy = 0, system = 4):
        file_name = os.path.basename(input_path).replace(".", "_")
        output_dir = "{}{}{}".format(output_dir, os.sep, file_name)
        try:
            os.makedirs(output_dir)
        except:
            pass
        
        self.output_path = output_dir + os.sep +file_name + ".dna"
        self.input_path = input_path
        self.split_len = split_len
        self.homo = homo
        self.gc_standard = gc_standard
        self.redanduncy = redanduncy
        self.system = system

    

    def transQuaternarySeq(self):
        self.quater_seq = []
        self.map_place = []




    
    def splitQuaternarySeq(self):


        # index 
        remainder = int(len(self.quater_seq)%self.split_len)
        if remainder == 0:
            seq_num = int(len(self.quater_seq)/self.split_len)
        else:
            seq_num = int(len(self.quater_seq)/self.split_len) +1
        idnex_len = len(tool().decimal2OtherSystem(seq_num-1,self.system))

        quater_seq_list = []
        index = 0
        for i in range(0, len(self.quater_seq), self.split_len):
            index_quater = tool().decimal2OtherSystem(index, self.system, idnex_len)

            # # transform index
            # index_quater = transSeqRandom(index_quater)

            seq = self.quater_seq[i: i+self.split_len]
            seq = index_quater + seq
            quater_seq_list.append(seq)
            index += 1

        # zero padding random
        if len(quater_seq_list[-1]) < self.split_len + idnex_len:
            padding_num = self.split_len+idnex_len-len(quater_seq_list[-1])

            ## add random int
            # set seed. I like this number, no reason!
            random.seed(3)
            for i in range(padding_num):
                quater_seq_list[-1] += str(random.randint(0, self.system-1))

        self.quater_seq_list = quater_seq_list
        self.idnex_len = idnex_len
        self.seq_num = len(quater_seq_list)


    def addRedundancy(self):
        if self.redanduncy == 0:
            return

        seq_num = len(self.quater_seq_list)*(self.redanduncy+1)
        if seq_num > int(seq_num):
            seq_num = int(seq_num) + 1
        else:
            seq_num = int(seq_num)

        # new index 
        idnex_len = len(tool().decimal2OtherSystem(seq_num-1,self.system))

        # add redanduncy
        redanduncy_list = []
        quater_seq_list_no_index = [i[self.idnex_len:] for i in self.quater_seq_list]
        step = int(1/self.redanduncy)

        for i in range(0, len(quater_seq_list_no_index), step):
            try:
                _seq_list = [quater_seq_list_no_index[j] for j in range(i, i+step)]
            except:
                _seq_list = [j for j in quater_seq_list_no_index[i:]]


            xor_seq = tool().xor(_seq_list)
            redanduncy_list.append(xor_seq)

        quater_seq_list_no_index += redanduncy_list

        # add index
        quater_seq_list = []
        for i in range(len(quater_seq_list_no_index)):
            index =  tool().decimal2OtherSystem(i, self.system, idnex_len)
            quater_seq = index + quater_seq_list_no_index[i]
            quater_seq_list.append(quater_seq)


        self.idnex_len = idnex_len
        self.quater_seq_list = quater_seq_list


    def homoControlWithSeqInsert(self):
        def splitSeq(seq):
            idnex_len = len(self.quater_seq_list[0]) - self.split_len
            part_seq_list = []
            index = seq[:idnex_len]
            info_seq = seq[idnex_len:]
            remainder = self.split_len%self.homo
            if remainder == 0:
                part_seq_len = int(self.split_len/self.homo)
            else:
                part_seq_len = int(self.split_len/self.homo)+1

            part_index_len = len(tool().decimal2OtherSystem(self.homo-1, self.system))
            for i in range(self.homo):
                part_index = tool().decimal2OtherSystem(i, self.system, part_index_len)

                part_seq = info_seq[i*part_seq_len: (i+1)*part_seq_len]
                part_seq += ["0"]*(part_seq_len-len(part_seq))
                part_seq = index + part_seq + part_index
                part_seq_list.append(part_seq)
            return part_seq_list

        def insertSeqAndJudge(part_seq, seq, step=self.homo):
            # insert base
            inserted_seq = []
            for i in range(len(part_seq)):
                inserted_seq += [part_seq[i]]
                inserted_seq += seq[i*step: (i+1)*step]
            inserted_seq += seq[(i+1)*step:] + [self.map_place[0]]

            # judge
            inserted_seq += ["$"]
            i, j = 0, 1
            while True:
                if i == len(inserted_seq)-1:
                    break
                if inserted_seq [j] == inserted_seq[i]:
                    j += 1
                else:
                    if j - i > self.homo:
                        return []
                    i = j 
                    j = i+1

            return inserted_seq[:-1]

        homo_control_list = []
        cannot_control_list = []
        quater_seq_list = copy.deepcopy(self.quater_seq_list)


        list_num_tag = len(quater_seq_list)
        loop_num = 0
        while True:
            # print(list_num_tag, loop_num)
            if list_num_tag == loop_num:
            # if len(quater_seq_list) < self.homo+1:
                # append quater seq to cannot control list
                for seq in quater_seq_list:
                    # "0123" at the start position is the mapping placeholder
                    # cannot_control_list.append(self.map_place + seq)
                    cannot_control_list.append(seq)
                break

            # list第一个元素进行切分，
            _homo_control_list = []
            quater_seq_found = []
            quater_seq_found_idnex = []

            part_seq_list = splitSeq(quater_seq_list[0])
            for part_seq in part_seq_list:
                for i in range(1, len(quater_seq_list)):

                    quater_seq = quater_seq_list[i]
                    # 已找到符合条件的序列，跳过
                    if i in quater_seq_found_idnex:
                        continue

                    inserted_seq = insertSeqAndJudge(part_seq, quater_seq)
                    if inserted_seq != []:
                        # quater_seq_found.append(part_seq)
                        quater_seq_found_idnex.append(i)
                        _homo_control_list.append(inserted_seq)
                        # _homo_control_list.append(inserted_seq+[self.map_place[0]])
                        break

                if len(_homo_control_list) == len(part_seq_list):
                    break

            if len(_homo_control_list) == len(part_seq_list):
                homo_control_list += _homo_control_list

                quater_seq_found_idnex.sort(reverse=True)
                for index in quater_seq_found_idnex:
                    del quater_seq_list[index]
                del quater_seq_list[0]

                list_num_tag = len(quater_seq_list)
                loop_num = 0
            else:
                # cannot_control_list.append(quater_seq_list[0])
                quater_seq_list = quater_seq_list[1:] + [quater_seq_list[0]]
                loop_num += 1

            # del quater_seq_list[0]
            
            
        self.homo_control_list = homo_control_list
        self.cannot_control_list = cannot_control_list

   

    def homoControlwithRandomInsert(self):
        if self.cannot_control_list == []:
            self.total_seq_list = self.homo_control_list
            return 

        try:
            add_num = len(self.homo_control_list[0]) - len(self.cannot_control_list[0])
        except:
            add_num = len(self.cannot_control_list[0])/self.homo
            if add_num > int(add_num):
                add_num = int(add_num) + 1
            else:
                add_num = int(add_num)
            # add tag position
            add_num += 1


        control_seq_list = []

        for quater_seq in self.cannot_control_list:
            control_seq = []
            for i in range(add_num):
                neibhor_base = []
                if control_seq != []:
                    neibhor_base.append(control_seq[-1])

                if quater_seq[i*self.homo: (i+1)*self.homo] != []:
                    neibhor_base.append(quater_seq[i*self.homo])

                ## insert
                # add random number
                for n in range(self.system):
                    if str(n) not in neibhor_base:
                        control_seq += [str(n)]
                        control_seq += quater_seq[i*self.homo: (i+1)*self.homo]
                        break

            # add 3 at the end of control_seq
            # add_num contains marker, minus 1
            control_seq = control_seq[:-1] + [self.map_place[-1]]
            control_seq_list.append(control_seq)

        self.total_seq_list = self.homo_control_list + control_seq_list

    

    def chooseMappingAndGetNtSeq(self):
        mapping_list = []
        if self.system > 4:
            return

        nt_seq_list = []
        for seq in self.total_seq_list:

            tag = seq[-1]
            # chose the suitable combination
            at_cha, cg_cha = [], []
            tag = len(seq)+1
            chose_combin = []
            num_dic = {str(i): seq.count(str(i)) for i in range(self.system)}
            standard_value = len(seq)*self.gc_standard

            for p in combinations(num_dic.values(), 2):
                _tag = abs(sum(p) - standard_value)
                if _tag < tag: 
                    tag = _tag
                    chose_combin = list(p)

            # confirm the map_dic
            for character, count in num_dic.items():
                if count in chose_combin:
                    cg_cha.append(character)
                    index = chose_combin.index(count)
                    chose_combin[index] = ""
                else:
                    at_cha.append(character)

            map_dic = {at_cha[0]: "A", at_cha[1]: "T", cg_cha[0]: "C", cg_cha[1]: "G"}
            key_list = [int(i) for i in map_dic.keys()]
            key_list.sort()
            _mapping = [map_dic[str(key)] for key in key_list]
            mapping_list.append("-".join(_mapping))

            # nt_seq = []
            # print(seq)
            # for num in seq:
            #     # print(num)
            #     nt_seq.append(map_dic[num])

            nt_seq = [map_dic[i] for i in seq]
            # change marker; uncontrol is AC, control is TG
            # if tag == "0":
            #     if nt_seq[-2] == "T":
            #         nt_seq[-1] = "G"
            #     else:
            #         nt_seq[-1] = "T"
            # elif tag == self.map_place[-1]:
            #     if nt_seq[-2] == "A":
            #         nt_seq[-1] = "C"
            #     else:
            #         nt_seq[-1] = "A"




            nt_seq_list.append(nt_seq)

        self.nt_seq_list = nt_seq_list
        self.mapping_list = mapping_list


    def addMappingID2Seq(self):
        if self.system != 4:
            return

        mapping_type_list = list(set(self.mapping_list))
        mapping_len = len(tool().decimal2OtherSystem(len(mapping_type_list), self.system))

        mapping_id_dic = {}
        for i in range(len(mapping_type_list)):
            trans_num = tool().decimal2OtherSystem(i, self.system, mapping_len)
            mapping = mapping_type_list[i]

            # 4-bases system
            # map_dic = {str(x): _mapping_l[x] for x in range(self.system)}
            map_dic = {"0":"A", "1":"C", "2":"G", "3":"T"}

            trans_num_seq = [map_dic[x] for x in trans_num]
            mapping_id_dic[mapping] = trans_num_seq

        # add mapping at start of seq
        for i in range(len(self.mapping_list)):
            mapping = self.mapping_list[i]
            seq = self.nt_seq_list[i]

            seq = mapping_id_dic[mapping] + seq 
            self.nt_seq_list[i] = seq 


        self.mapping_type_list = mapping_type_list




    def saveNtSeq(self):
        if not self.output_path.endswith(".dna"):
            self.output_path += ".dna"

        if self.system ==  4:
            with open(self.output_path, "w") as f:
                for seq in self.nt_seq_list:
                    f.write("".join(seq)+"\n")

        else:
            with open(self.output_path, "w", encoding = "utf-8") as f:
                for seq in self.total_seq_list:
                    f.write(",".join(seq)+"\n")


    def saveMapping(self):
        if self.system != 4:
            return

        mapping_path = self.output_path+".mapping.dna"
        mapping_seq = "".join(self.mapping_type_list)
        mapping_seq = mapping_seq.replace("-", "")

        with open(mapping_path, "w") as f:
            f.write(mapping_seq+"\n")



    def saveConfig(self):
        
        config_path = self.output_path + ".config"
        with open(config_path, "w") as f:
            outline = "storageSystem\t{}\nfileExtension\t{}\nsplitLength\t{}\nrunningLength\t{}\nindexLength\t{}\nsequenceNumber\t{}\nredanduncy\t{}\n".format(
                      self.system, self.input_path.split(".")[-1], self.split_len, self.homo, self.idnex_len, self.seq_num, self.redanduncy)
            f.write(outline)

    def saveConfigIdentity(self):
        pass


    def main(self):
        self.transQuaternarySeq()
        self.splitQuaternarySeq()
        self.addRedundancy()
        self.homoControlWithSeqInsert()
        self.homoControlwithRandomInsert()
        self.chooseMappingAndGetNtSeq()
        self.addMappingID2Seq()
        self.saveNtSeq()
        self.saveMapping()
        self.saveConfig()
        self.saveConfigIdentity()

class textEncode(insertEncode):

    def transQuaternarySeq(self):

        map_place = [str(i) for i in range(self.system)]
        ch_num_list = []
        with open(self.input_path, encoding = "utf-8") as f:
            for line in f:
                for ch in line:
                    ch_num = ord(ch)
                    ch_num_list.append(ch_num)

        max_num = max(ch_num_list)
        bits = math.log(max_num, 2)
        if int(bits) < bits:
            bits = int(bits) + 1


        file_system = 2**bits


        
        file_bits_num, unit_len = tool().getLCM(file_system, self.system)
        # print(file_bits_num, unit_len)

        # get quater_seq
        quater_seq = []
        for i in range(0, len(ch_num_list), file_bits_num):
            num_unit = ch_num_list[i: i+file_bits_num]
            num_unit += [0]*(file_bits_num-len(num_unit))
            decimal_number = tool().convert2Decimal(num_unit, file_system)
            quater_ch_num = tool().decimal2OtherSystem(decimal_number, self.system, unit_len)
            quater_seq += quater_ch_num

        self.map_place = map_place
        self.quater_seq = quater_seq
        self.file_system = file_system
        self.ch_num = len(ch_num_list)



    def saveConfigIdentity(self):
        config_path = self.output_path+ ".config"
        density = os.path.getsize(self.input_path)*8/(len(self.total_seq_list[0])*len(self.total_seq_list))
        
        with open(config_path, "a") as f:
            outline = "density\t{}\nfileSystem\t{}\ncharacterNumber\t{}\n".format(density, self.file_system, self.ch_num)
            f.write(outline)


class textEncode_unicode(insertEncode):
    # data file is coded by utf-8

 

    def transQuaternarySeq(self):

        def convertQuaternaryNum(ch_num):
            pre_key= 0
            bytes_num_dic, cha = tool().getTxtCodingDic(self.system)
            for key, value in bytes_num_dic.items():
                if ch_num >= key:
                    pre_key = key
                    continue

                _ch_num = ch_num - pre_key
                if self.system < 256:
                    num_list = tool().decimal2OtherSystem(_ch_num, self.system, precision=value*cha)
                else:
                    num_list = tool().decimal2OtherSystem(_ch_num, self.system, precision=2*value-1)

                # marker change
                for i in range(value-1):
                    num_list[i] = str(self.system-1)
                break

            return num_list


        map_place = [str(i) for i in range(self.system)]

        quater_seq = []
        with open(self.input_path, encoding="utf-8") as f:
            for line in f:
                for ch in line:
                    ch_num = ord(ch)
                    quater_ch_num = convertQuaternaryNum(ch_num)
                    quater_seq += quater_ch_num

        self.quater_seq = quater_seq
        self.map_place = map_place

    def saveConfigIdentity(self):
        config_path = self.output_path+ ".config"
        density = os.path.getsize(self.input_path)*8/(len(self.total_seq_list[0])*len(self.total_seq_list))
        
        with open(config_path, "a") as f:
            outline = "density\t{}\n".format(density)
            f.write(outline)



class imageEncode(insertEncode):
    def transQuaternarySeq(self):
        quater_seq = []
        # CLM
        file_bits_num, precision = tool().getLCM(256, self.system)

        im = Image.open(self.input_path)
        try:
            dpi = list(im.info.get("dpi"))[0]
        except:
            dpi = 200

        pix = im.load()
        width, height = im.size
        rgb_list = []
        for y in range(height):
            for x in range(width):
                r, g, b = pix[x, y]
                rgb_list += [r,g,b]

        # get quater_seq
        for i in range(0, len(rgb_list), file_bits_num):
            num_unit = rgb_list[i: i+file_bits_num]
            num_unit += ["0"]*(file_bits_num-len(num_unit))
            decimal_number = tool().convert2Decimal(num_unit, 256)
            quater_seq += tool().decimal2OtherSystem(decimal_number, self.system, precision)


        map_place = [str(i) for i in range(self.system)]

        self.quater_seq = quater_seq
        self.width = width
        self.height = height
        self.dpi = dpi
        self.map_place = map_place

    def saveConfigIdentity(self):
        config_path = self.output_path+ ".config"
        density = self.width*self.height*3*8/(len(self.total_seq_list[0])*len(self.total_seq_list))

        with open(config_path, "a") as f:
            outline = "width\t{}\nheight\t{}\ndpi\t{}\ndensity\t{}\n".format(self.width, self.height, self.dpi, density)
            f.write(outline)


class audioEncode(insertEncode):
    
    def transQuaternarySeq(self):
        sound = AudioSegment.from_file(self.input_path)
        samples = sound.get_array_of_samples()
        shifted_sample = list(samples)

        add_num = int(sound.max_possible_amplitude)
        shifted_sample = [i+add_num for i in shifted_sample]

        # unit_len
        file_system = int(2*sound.max_possible_amplitude)
        file_bits_num, unit_len = tool().getLCM(file_system, self.system)

        quater_seq = []
        for i in range(0, len(shifted_sample), file_bits_num):
            sample_unit = shifted_sample[i: i+file_bits_num]
            sample_unit += [0]*(file_bits_num-len(sample_unit))

            decimal_number = tool().convert2Decimal(sample_unit, file_system)
            quater_seq += tool().decimal2OtherSystem(decimal_number, self.system, unit_len)


        map_place = [str(i) for i in range(self.system)]


        self.quater_seq = quater_seq
        self.max = int(sound.max_possible_amplitude)
        self.sample_num = len(shifted_sample)
        self.channels = sound.channels
        self.array_type = sound.array_type
        self.frame_rate = sound.frame_rate
        self.shifted_sample = shifted_sample

        self.map_place = map_place

    def saveConfigIdentity(self):
        config_path = self.output_path+ ".config"
        density = math.log(self.max, 2)*self.sample_num/(len(self.total_seq_list[0])*len(self.total_seq_list))

        with open(config_path, "a") as f:
            outline = "maxPossibleAmplitude\t{}\nsampleNumber\t{}\nchannels\t{}\narrayType\t{}\nframeRate\t{}\ndensity\t{}\n".format(
                        self.max, self.sample_num, self.channels, self.array_type, self.frame_rate, density)
            f.write(outline)


class videoEncode:
    def __init__(self, input_path, output_dir, split_len=200, homo=4, gc_standard = 0.5, redanduncy = 0, system = 4):
        file_name = os.path.basename(input_path).replace(".", "_")
        output_dir = "{}{}{}".format(output_dir, os.sep, file_name)
        tmp_dir = "{}{}{}_encodeTmp".format(output_dir, os.sep, file_name)
        try:
            os.makedirs(tmp_dir)
        except:
            pass

        self.input_path = input_path
        self.split_len = split_len
        self.homo = homo
        self.gc_standard = gc_standard
        self.tmp_dir = tmp_dir
        self.output_dir = output_dir
        self.file_name = file_name
        self.redanduncy = redanduncy
        self.system = system


    def getAudio(self):
        audio_path = "{}{}{}.mp3".format(self.tmp_dir, os.sep, self.file_name)
        clip = mp.VideoFileClip(self.input_path)
        audio = clip.audio

        # slience video
        if audio == None:
            self.audio_path = None
            return

        audio.write_audiofile(audio_path)
        clip.close()

        self.audio_path = audio_path


    def getImages(self):
        image_path_list = []
        cap = cv2.VideoCapture(self.input_path)
        fps = cap.get(cv2.CAP_PROP_FPS)


        count = 0    
        while True:
            flag, frame = cap.read()
            if flag == False:
                break

            # no Chinese
            image_path = "{}{}{}_{}.jpg".format(self.tmp_dir, os.sep, self.file_name, count)
            image_path_list.append(image_path)
            cv2.imwrite(image_path, frame)

            count += 1

        cap.release()



        self.image_path_list = image_path_list
        self.fps = int(fps)


    def codingFiles(self):

        if self.audio_path == None:
            input_path_list = self.image_path_list
        else:
            input_path_list = [self.audio_path] + self.image_path_list

        params = [(input_path, self.output_dir, self.split_len, self.homo, self.gc_standard, self.redanduncy, self.system) for input_path in input_path_list]

        # default 16 processing
        with Pool(processes = 16) as pool:
            results = pool.starmap(codingFiles_worker, params)

       
        # density
        bits_total, bases_total = 0, 0
        for bits, bases in results:
            bits_total += bits
            bases_total += bases
        density = bits_total/bases_total

        self.density = density





    def saveConfig(self):
        # video config
        config_path = self.output_dir + os.sep + self.file_name + ".config"
        outline = "fileExtension\t{}\nFPS\t{}\ndensity\t{}\n".format(self.input_path.split(".")[-1], self.fps, self.density)
        with open(config_path, "w") as f:
            f.write(outline)


    def main(self):
        self.getAudio()
        self.getImages()
        self.codingFiles()
        self.saveConfig()


## video encoding multiprocessing function, global needed
def codingFiles_worker(input_path, output_dir, split_len, homo, gc_standard, redanduncy, system):

    print('Worker (PID):', os.getpid(), ', File:', input_path)

    if input_path.endswith(".mp3"):
        au_cl = audioEncode(input_path, output_dir, split_len, homo, gc_standard, redanduncy, system)
        au_cl.main()

        bits = math.log(au_cl.max, 2) * au_cl.sample_num
        bases = len(au_cl.total_seq_list[0]) * len(au_cl.total_seq_list)
    
    else:
        im_cl = imageEncode(input_path, output_dir, split_len, homo, gc_standard, redanduncy, system)
        im_cl.main()

        bits = im_cl.width * im_cl.height * 3 * 8
        bases = len(im_cl.total_seq_list[0]) * len(im_cl.total_seq_list)

    return [bits, bases]




class insertDecode:
    def __init__(self, seq_dir, save_dir):
        try:
            os.makedirs(save_dir)
        except:
            pass 


        file_list = os.listdir(seq_dir)
        for file in file_list:
            if file.endswith(".dna"):
                if not "mapping" in file:
                    seq_path = "{}{}{}".format(seq_dir, os.sep, file)
                    save_path = "{}{}{}".format(save_dir, os.sep, file)
                    break
            else:
                seq_path = ""
                save_path = ""

        self.seq_path = seq_path
        self.save_path = save_path


    def readConfig(self):
        config_path = self.seq_path + ".config"
        config_dic = {}

        with open(config_path) as f:
            for i in f:
                key, value = i.strip().split("\t")
                config_dic[key] = value

        # correct the file extension
        if not self.save_path.endswith(config_dic["fileExtension"]):
            self.save_path += ".{}".format(config_dic["fileExtension"])

        self.config_dic = config_dic
        self.system = int(self.config_dic["storageSystem"])

    def readSeq(self):
        nt_seq_list = []
        if self.system == 4:
            with open(self.seq_path) as f:
                for seq in f:
                    nt_seq_list.append(list(seq.strip()))

        else:
            with open(self.seq_path) as f:
                for seq in f:
                    nt_seq_list.append(seq.strip().split(","))

        self.nt_seq_list = nt_seq_list


    def readMappingSeq(self):
        if self.system != 4:
            return 

        mapping_path = self.seq_path + ".mapping.dna"
        with open(mapping_path) as f:
            seq = f.read().strip()

        map_list = []
        for i in range(0, len(seq), self.system):
            mapping = seq[i: i+self.system]
            mapping_dic = {mapping[x]:str(x) for x in range(self.system)}
            map_list.append(mapping_dic)

        self.map_list = map_list



    def recoverNtSeq(self):

        def transIndex2Int(num_seq):
                # int index
                for i in range(int(self.config_dic["indexLength"])):
                    num_seq[i] = int(num_seq[i])
                return num_seq


        remove_insert_len = int(self.config_dic["indexLength"]) + int(self.config_dic["splitLength"])
        # has 1 tag at the end of seq
        isnert_len = len(self.nt_seq_list[0])-1 - remove_insert_len

        quater_seq_list, part_seq_list = [], []
        if self.system == 4:
            map_len = len(tool().decimal2OtherSystem(len(self.map_list), self.system))
            map_trans_dic = {"A":"0", "C":"1", "G":"2", "T":"3"}

        for seq in self.nt_seq_list:
            original_seq = []
            insert_seq = []
            tag = seq[-1]

            if self.system == 4:
                map_id = int("".join([map_trans_dic[x] for x in seq[:map_len]]))
                nt_dic = self.map_list[map_id]
                seq = seq[map_len:-1]

            
            else:
                seq = seq[:-1]

            # if (tag == "0") or (tag in ["T", "G"]):
            #     isnert_len = control_isnert_len
            #     remove_insert_len = control_remove_insert_len

            # else:
            #     isnert_len = uncontrol_isnert_len
            #     remove_insert_len = uncontrol_remove_insert_len




            step = int(self.config_dic["runningLength"])+1
            for i in range(0, len(seq), step):
                unit_seq = seq[i: i+step]

                if len(insert_seq) == isnert_len:
                    original_seq += unit_seq
                    continue
                if len(original_seq) == remove_insert_len:
                    insert_seq += unit_seq
                    continue

                insert_seq += [unit_seq[0]]
                for nt in unit_seq[1:]:
                    if len(original_seq) == remove_insert_len:
                        insert_seq += [nt]
                    else:
                        original_seq += [nt] 


            # trans to quaternary seq
            if self.system == 4:
                quater_seq = [nt_dic[i] for i in original_seq]
                quater_seq = transIndex2Int(quater_seq)
                quater_seq_list.append(quater_seq)

                if nt_dic[tag] == "0":
                    part_seq = [nt_dic[i] for i in insert_seq]
                    part_seq_list.append(part_seq)

            else:
                quater_seq = transIndex2Int(original_seq)
                quater_seq_list.append( quater_seq)


                if tag == "0":
                    part_seq_list.append(insert_seq)



            # if (tag == "0") or (tag in ["T", "G"]):
            #     nt_dic = {original_seq[i]: str(i) for i in range(self.system)}
            #     quater_seq = [nt_dic[i] for i in original_seq][self.system:]  
            #     quater_seq = transIndex2Int(quater_seq)
            #     quater_seq_list.append(quater_seq)

            # else:
            #     if self.system == 4:
            #         nt_dic = {"A": "0", "T": "1", "C": "2", "G": "3"}
            #         quater_seq = [nt_dic[i] for i in original_seq]
            #         quater_seq_list.append(quater_seq)
            #         quater_seq = transIndex2Int(quater_seq)
            #     else:
            #         original_seq = transIndex2Int(original_seq)
            #         quater_seq_list.append(original_seq)

            # if (tag == "0") or (tag in ["T", "G"]):
            #     part_seq = [nt_dic[i] for i in insert_seq]
            #     part_seq_list.append(part_seq)
                # quater_seq_part_seq_dic[quater_seq] = part_seq
            # else:
            #     quater_seq_part_seq_dic[quater_seq] = ""


        # self.quater_seq_part_seq_dic = quater_seq_part_seq_dic
        quater_seq_list.sort()
        part_seq_list.sort()
        self.quater_seq_list = quater_seq_list
        self.part_seq_list = part_seq_list



    def recoverPartSeq(self):
        part_recover_seq_list = []    
        part_index_len = len(tool().decimal2OtherSystem(int(self.config_dic["runningLength"])-1, self.system))

        for i in range(0, len(self.part_seq_list), int(self.config_dic["runningLength"])):
            index = [int(x) for x in self.part_seq_list[i][:int(self.config_dic["indexLength"])]]

            unit_seq_list = self.part_seq_list[i: i+ int(self.config_dic["runningLength"])]
            index_head_unit_seq_list = [seq[-1*part_index_len:] + seq[len(index):-1*part_index_len] for seq in unit_seq_list]
            index_head_unit_seq_list.sort()

            part_recover_seq = []
            for seq in index_head_unit_seq_list:
                part_recover_seq += seq[part_index_len:]

            part_recover_seq = index + part_recover_seq
            part_recover_seq_list.append(part_recover_seq)

        self.part_recover_seq_list = part_recover_seq_list
        self.total_seq_list = part_recover_seq_list + self.quater_seq_list
        self.total_seq_list.sort()


    def redanduncyCorrection(self):
        # 需要根据实际的测序数据再详细写，目前直接截断了冗余部分
        if self.config_dic["redanduncy"] == "0":
            return

        total_seq_list = self.total_seq_list
        total_seq_list = total_seq_list[:int(self.config_dic["sequenceNumber"])]

        self.total_seq_list = total_seq_list

    def combineSeq(self):
        total_seq_list = self.total_seq_list
        quater_seq_complement = []

        for i in range(len(total_seq_list)):
            quater_seq = total_seq_list[i]

            # remove zero padding
            quater_seq = quater_seq[: int(self.config_dic["indexLength"]) + int(self.config_dic["splitLength"])]
            # with index
            total_seq_list[i] = quater_seq

            quater_seq = quater_seq[int(self.config_dic["indexLength"]):]
            quater_seq_complement += quater_seq

        self.total_seq_list = total_seq_list
        self.quater_seq_complement = quater_seq_complement



    def recoverFile(self):
        pass



    def writeFile(self):
        pass



    def main(self):
        self.readConfig()
        self.readSeq()
        self.readMappingSeq()
        self.recoverNtSeq()
        self.recoverPartSeq()
        self.redanduncyCorrection()
        self.combineSeq()
        self.recoverFile()
        self.writeFile()


class textDecode(insertDecode):


    def recoverFile(self):
        file_system = int(self.config_dic["fileSystem"])
        file_bits_num, unit_len = tool().getLCM(file_system, self.system)

        # remove padding
        character_num = int(np.ceil(int(self.config_dic["characterNumber"])/file_bits_num)*unit_len)
        self.quater_seq_complement = self.quater_seq_complement[:character_num]

        cha_list = []
        for i in range(0, len(self.quater_seq_complement), unit_len):
            quater_ch_num = self.quater_seq_complement[i: i+unit_len]
            decimal_number = tool().convert2Decimal(quater_ch_num, self.system)
            cha_num = tool().decimal2OtherSystem(decimal_number, file_system, file_bits_num)
            _cha_list = [chr(int(x)) for x in cha_num]
            cha_list += _cha_list


        self.cha_seq = "".join(cha_list)


    def writeFile(self):
        with open(self.save_path, "w", encoding="utf-8") as f:
            f.write(self.cha_seq)


class textDecode_unicode(insertDecode):

    def removeTotalZeroPadding(self):
        while True:
            if self.quater_seq_complement[-1] == "0":
                self.quater_seq_complement = self.quater_seq_complement[:-1]
            else:
                break


    def recoverFile(self):
        bytes_num_dic, cha = tool().getTxtCodingDic(self.system)
        bytes_num_dic = {y: x for x,y in bytes_num_dic.items()}
        bytes_num_dic[0] = 0
        # bytes_num_dic = {0:0, 1:192, 2:12480, 3:798912, 4:51130560}

        # 连续最高位的个数决定了character的个数，4进制中最高位为3
        cha_list = []
        p1 = 0
        while True:
            if p1 > len(self.quater_seq_complement)-1:
                break

            if int(self.quater_seq_complement[p1]) < self.system-1:
                unit_num = 1
                p2 = p1

            else:
                # find consecutive maximum number
                p2 = p1+1
                while True:
                    if self.quater_seq_complement[p2] != self.quater_seq_complement[p1]:
                        break
                    p2 += 1
                unit_num = p2 - p1 +1

            if self.system < 256:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                cha_num = self.quater_seq_complement[p2: p1+cha*unit_num]
                cha_num = tool().convert2Decimal(cha_num, self.system) + bytes_num_dic[unit_num-1]
                cha_list.append(chr(cha_num))
                p1 += cha*unit_num

            else:
                cha_num = self.quater_seq_complement[p2: p2+cha*unit_num]
                cha_num = tool().convert2Decimal(cha_num, self.system) + bytes_num_dic[unit_num-1]
                cha_list.append(chr(cha_num))
                p1 += 2*unit_num-1



        self.cha_seq = "".join(cha_list)

    def writeFile(self):
        with open(self.save_path, "w", encoding="utf-8") as f:
            f.write(self.cha_seq)


    
class imageDecode(insertDecode):

    def recoverFile(self):
        file_bits_num, precision = tool().getLCM(256, self.system)
        rgb_list = []
        for i in range(0, len(self.quater_seq_complement), precision):
            num = self.quater_seq_complement[i: i+precision]
            decimal_number = tool().convert2Decimal(num, self.system)
            rgb_num = [int(x) for x in tool().decimal2OtherSystem(decimal_number, 256, file_bits_num)]
            rgb_list += rgb_num
        self.rgb_list = rgb_list


    def writeFile(self):
        width = int(self.config_dic["width"])
        height = int(self.config_dic["height"])
        dpi = int(self.config_dic["dpi"])

        im = Image.new("RGB", size= (width, height))


        n = 0
        for _h in range(height):
            for _w in range(width):
                _color = tuple(self.rgb_list[n: n+3])
                _im = Image.new("RGB", size=(1,1), color = _color)
                im.paste(_im, (_w, _h, _w+1, _h+1))
                n+=3
        im.save(self.save_path, dpi = (dpi, dpi))
        im.close()

    # def md5Checkout(self):
    #     # 恢复的图片压缩率可能不一致 但是不影响信息传递  此处不进行校验
    #     pass


class audioDecode(insertDecode):
    def recoverFile(self):
        max_possible_amplitude = int(self.config_dic["maxPossibleAmplitude"])
        sample_num = int(self.config_dic["sampleNumber"])

        # get unit_len
        file_system = int(max_possible_amplitude*2)
        file_bits_num, unit_len = tool().getLCM(file_system, self.system)

        
        sample_list = []
        for i in range(0, sample_num, file_bits_num):
            unit = self.quater_seq_complement[i*unit_len: (i+1)*unit_len]
            decimal_number = tool().convert2Decimal(unit, self.system)
            unit_num = [int(x) for x in tool().decimal2OtherSystem(decimal_number, file_system, file_bits_num)]
            sample_list+= unit_num
        sample_list = sample_list[:sample_num]

        # to array
        dtype = "int{}".format(int(math.log(max_possible_amplitude*2,2)))
        sample_array = np.array(sample_list, dtype = dtype)
        sample_array -= int(max_possible_amplitude)
        sample_array = array.array(self.config_dic["arrayType"], sample_array)

        self.sample_array = sample_array
        self.dtype = dtype


    def writeFile(self):
        # to AudioSegment
        empty_array = np.array([], dtype=self.dtype)

        wav_io = io.BytesIO()
        scipy.io.wavfile.write(wav_io, int(self.config_dic["frameRate"]), empty_array)
        wav_io.seek(0)
        empty_sound = AudioSegment.from_wav(wav_io)

        empty_sound = empty_sound.set_channels(int(self.config_dic["channels"]))
        empty_sound = empty_sound.set_frame_rate(int(self.config_dic["frameRate"]))
        sound = empty_sound._spawn(self.sample_array)

        sound.export(self.save_path, format(self.config_dic["fileExtension"]))

        # self.sound = sound




class videoDecode:

    def __init__(self, seq_dir, save_dir):
        tmp_dir = "{}{}{}_decodeTmp".format(save_dir, os.sep, os.path.basename(seq_dir))
        try:
            os.makedirs(tmp_dir)
        except:
            pass 


        self.seq_dir = seq_dir
        self.save_path = save_dir + os.sep + os.path.basename(seq_dir) + ".mp4"
        self.tmp_dir = tmp_dir

    
    def recoverFiles(self):
        files_dir = ["{}{}{}".format(self.seq_dir, os.sep, i) for i in os.listdir(self.seq_dir) if (i.endswith("_mp3") or i.endswith("_jpg"))]
        params = [(i, self.tmp_dir) for i in files_dir]

        with Pool(processes = 16) as pool:
            results = pool.starmap(recoverFiles_worker, params)

        
        for i in results:
            if i.endswith(".mp3"):
                self.audio_path = i
                self.audio_judge = True
                break
            
            else:
                self.audio_judge = False


    def recoverVideo(self):
        # get video config
        config_path = self.seq_dir + os.sep + os.path.basename(self.seq_dir) + ".config"
        config_dic = {}
        with open(config_path) as f:
            for line in f:
                key, value = line.strip().split("\t")
                config_dic[key] = value

        # recover images video
        fps = int(config_dic["FPS"])
        iamge_list = [i for i in os.listdir(self.tmp_dir) if i.endswith("jpg")]
        with imageio.get_writer(self.save_path, fps=fps) as writer:


            image_num = len(os.listdir(self.tmp_dir)) -1
            for i in range(image_num):
                for iamge in iamge_list:
                    if "_{}_".format(i) in iamge:
                        image_path = "{}/{}".format(self.tmp_dir, iamge)
                        break
                img = imageio.imread(image_path)
                writer.append_data(img)

    def combineVideoAudio(self):
        # slience video
        if self.audio_judge == False:
            return 

        video = mp.VideoFileClip(self.save_path)
        audio = mp.AudioFileClip(self.audio_path)
        video = video.set_audio(audio)

            
        video.write_videofile(self.save_path)




    def main(self):
        self.recoverFiles()
        self.recoverVideo()
        self.combineVideoAudio()


def recoverFiles_worker(seq_dir, file_dir):

    print('Worker (PID):', os.getpid(), ', File:', seq_dir)

    if seq_dir.endswith("_mp3"):
        audioDecode(seq_dir, file_dir).main()
        audio_path = ["{}/{}".format(file_dir, i) for i in os.listdir(file_dir) if i.endswith("mp3")][0]
        return audio_path

    else:
        imageDecode(seq_dir, file_dir).main()
        return "0"



###################  main function #############################

def encode(input_path, encode_dir, split_len=200, homo=4, gc_standard = 0.5, redanduncy = 0, system = 4):
    file_extension = os.path.basename(input_path).split(".")[-1]

    if file_extension in ["txt"]:
        cl = textEncode(input_path, encode_dir, split_len = split_len, homo = homo, gc_standard = gc_standard, redanduncy = redanduncy, system = system)
    elif file_extension in ["jpg", "jpeg", "png"]:
        cl = imageEncode(input_path, encode_dir, split_len = split_len, homo = homo, gc_standard = gc_standard, redanduncy = redanduncy, system = system)
    elif file_extension in ["mp3", "wav"]:
        cl = audioEncode(input_path, encode_dir, split_len = split_len, homo = homo, gc_standard = gc_standard, redanduncy = redanduncy, system = system)
    elif file_extension in ["mp4"]:
        cl = videoEncode(input_path, encode_dir, split_len = split_len, homo = homo, gc_standard = gc_standard, redanduncy = redanduncy, system = system)
    else:
        raise ValueError("The type of {} file can't be processed!".format(file_extension))

    cl.main()

    return cl


def decode(encode_dir, decode_dir):
    # read config_dic 
    config_path = ["{}/{}".format(encode_dir, i) for i in os.listdir(encode_dir) if i.endswith("config")][0]
    with open(config_path) as f:
        for i in f:
            if i.startswith("fileExtension"):
                file_extension = i.strip().split("\t")[1]
                break


    if file_extension in ["txt"]:
        cl = textDecode(encode_dir, decode_dir)
    elif file_extension in ["jpg", "jpeg", "png"]:
        cl = imageDecode(encode_dir, decode_dir)
    elif file_extension in ["mp3", "wav"]:
        cl = audioDecode(encode_dir, decode_dir)
    elif file_extension in ["mp4"]:
        cl = videoDecode(encode_dir, decode_dir)
    else:
        raise ValueError("The config file is not correct!\nFile extension: {}".format(file_extension))

    cl.main()

    return cl



def linuxComand():
    parser = argparse.ArgumentParser(description='Direct codec algorithm of DNA data storage')

    subparsers = parser.add_subparsers(dest='command')
    
    add_parser = subparsers.add_parser('encode')
    add_parser.add_argument('-i', '--input_file', type=str, required=True, help='Input data file')
    add_parser.add_argument('-o_d', '--output_dir', type=str, required=True, help='Result directory')
    add_parser.add_argument('-l', '--split_len', type=int, default = 200, help='The length of original segment, not the final length, default=[200]')
    add_parser.add_argument('-homo', '--homopolymer', type=int, default = 4, help='The length of homopolymers, default = [4]')
    add_parser.add_argument('-gc', '--gc_standard', type=float, default = 0.5, help='The expected standard GC ratio, default = [0.5]')
    add_parser.add_argument('-r', '--redanduncy', type=float, default = 0, help='Default = [0]')
    add_parser.add_argument('-s', '--system', type=int, default = 4, help='The type of bases of storage system, default = [4]')

    add_parser = subparsers.add_parser('decode')
    add_parser.add_argument('-i_d', '--input_dir', type=str, required=True, help='Input data directory')
    add_parser.add_argument('-o_d', '--output_dir', type=str, required=True, help='Result directory')

    args = parser.parse_args()


    # chose mode to run
    if args.command == 'encode':
        encode(args.input_file, args.output_dir, split_len=args.split_len, homo=args.homopolymer, gc_standard = args.gc_standard, 
                redanduncy = args.redanduncy, system = args.system)
    elif args.command == 'decode':
        decode(args.input_dir, args.output_dir)
    else:
        parser.print_help()  # 如果没有指定子命令，则打印帮助信息


if __name__ == "__main__":

    linuxComand()

    # input_path = "D:/BaiduSyncdisk/中科院先进院合成所/项目/2021_12_06-分割插入编码/修饰碱基测序方法构建/Nanopore/input/DNA.txt"
    # encode_dir = "D:/BaiduSyncdisk/中科院先进院合成所/项目/2021_12_06-分割插入编码/修饰碱基测序方法构建/Nanopore/encode"
    # # input_path = "D:/BaiduSyncdisk/中科院先进院合成所/项目/2021_12_06-分割插入编码/test_data/tt.txt"
    # # input_path = "D:/BaiduSyncdisk/中科院先进院合成所/项目/2021_12_06-分割插入编码/test_data/audio_test/test_1000.mp3"
    # # encode_dir = os.path.dirname(input_path)
    # redanduncy = 0
    # system = 4
    # split_len = 370
    # homo = 3
    # cl = encode(input_path, encode_dir, redanduncy = redanduncy, system = system, split_len = split_len, homo = homo)

    # quater_seq=cl.quater_seq
    # quater_seq_list = cl.quater_seq_list
    # homo_control_list = cl.homo_control_list
    # cannot_control_list = cl.cannot_control_list
    # total_seq_list = cl.total_seq_list
    # nt_seq_list = cl.nt_seq_list
    # shifted_sample = cl.shifted_sample

    # # gc_list = [(i.count("C")+i.count("G"))/len(i) for i in nt_seq_list]
    # # gc_list.sort()
    # # print(round(gc_list[0], 4), round(gc_list[-1],4))

    # #########################  decode
    # encode_dir = os.path.dirname(input_path)+os.sep+os.path.basename(input_path).replace(".", "_")
    # decode_dir = encode_dir
    # encode_dir += os.sep + "DNA_txt"
    # decode_dir = "D:/BaiduSyncdisk/中科院先进院合成所/项目/2021_12_06-分割插入编码/修饰碱基测序方法构建/Nanopore/decode"

    # cld = decode(encode_dir, decode_dir)


    # quater_seq_list_d = cld.quater_seq_list 
    # part_seq_list = cld.part_seq_list
    # part_recover_seq_list = cld.part_recover_seq_list
    # total_seq_list_d = cld.total_seq_list
    # quater_seq_complement = cld.quater_seq_complement
    # # nt_seq_list_d = cld.nt_seq_list
    # # sound = cld.sound



