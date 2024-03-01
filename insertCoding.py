# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 11:16:45 2024

@author: Wei Qiang
"""
import sys, copy, math, array, io, os, cv2
from itertools import combinations
from md5hash import scan
from PIL import Image
from pydub import AudioSegment
import numpy as np
import scipy.io.wavfile
import moviepy.editor as mp
import imageio

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


class insertEncode:

    def __init__(self, input_path, output_dir, split_len=200, homo=4, gc_standard = 0.5):
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

    

    def transQuaternarySeq(self):
        self.quater_seq = ""


    def randomQuaternarySeq(self):
        new_quater_seq = []

        tag = False
        s, e = 0, len(self.quater_seq)-1
        while True:
            if tag == True:
                break
            if s == e:
                new_quater_seq.append(self.quater_seq[s])
                tag = True
            else:
                new_quater_seq.append(self.quater_seq[s])
                new_quater_seq.append(self.quater_seq[e])
                if e == s+1:
                    tag = True
            s+= 1
            e-= 1

        self.quater_seq = "".join(new_quater_seq)

    def randomQuaternarySeq2(self):
        new_quater_seq = []
        for i in range(len(self.quater_seq)):
            num = int(self.quater_seq[i])
            add_num = i%4
            sum_r = (add_num + num)%4
            new_quater_seq.append(str(sum_r))

        self.quater_seq = "".join(new_quater_seq)


    
    def splitQuaternarySeq(self):

        # def transSeqRandom(num_seq):
        #     new_num_seq = []
        #     for i in range(len(num_seq)):
        #         num = int(num_seq[i])
        #         add_num = i%4
        #         sum_r = (add_num + num)%4
        #         new_num_seq.append(str(sum_r))

        #     new_num_seq = "".join(new_num_seq)
        #     return new_num_seq

        # index 
        remainder = int(len(self.quater_seq)%self.split_len)
        if remainder == 0:
            seq_num = int(len(self.quater_seq)/self.split_len)
        else:
            seq_num = int(len(self.quater_seq)/self.split_len) +1
        idnex_len = len(tool().decimal2OtherSystem(seq_num-1,4))

        quater_seq_list = []
        index = 0
        for i in range(0, len(self.quater_seq), self.split_len):
            index_quater = "".join(tool().decimal2OtherSystem(index, 4, idnex_len))

            # # transform index
            # index_quater = transSeqRandom(index_quater)

            seq = self.quater_seq[i: i+self.split_len]
            seq = index_quater + seq
            quater_seq_list.append(seq)
            index += 1

        # zero padding
        if len(quater_seq_list[-1]) < self.split_len + idnex_len:
            quater_seq_list[-1] += (self.split_len+idnex_len-len(quater_seq_list[-1]))*"0"

        self.quater_seq_list = quater_seq_list
        self.idnex_len = idnex_len







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

            part_index_len = len(tool().decimal2OtherSystem(self.homo-1, 4))
            for i in range(self.homo):
                part_index = "".join(tool().decimal2OtherSystem(i, 4, part_index_len))

                part_seq = info_seq[i*part_seq_len: (i+1)*part_seq_len]
                part_seq += "0"*(part_seq_len-len(part_seq))
                part_seq = index + part_seq + part_index
                part_seq_list.append(part_seq)
            return part_seq_list

        def insertSeqAndJudge(part_seq, seq, step=self.homo):
            # insert base
            inserted_seq = ""
            for i in range(len(part_seq)):
                inserted_seq += part_seq[i]
                inserted_seq += seq[i*step: (i+1)*step]

            # judge
            inserted_seq += "$"
            i, j = 0, 1
            while True:
                if i == len(inserted_seq)-1:
                    break
                if inserted_seq [j] == inserted_seq[i]:
                    j += 1
                else:
                    if j - i > self.homo:
                        return ""
                    i = j 
                    j = i+1

            return inserted_seq[:-1]

        homo_control_list = []
        cannot_control_list = []
        quater_seq_list = copy.deepcopy(self.quater_seq_list)


        list_num_tag = len(quater_seq_list)
        loop_num = 0
        while True:
            print(list_num_tag, loop_num)
            if list_num_tag == loop_num:
            # if len(quater_seq_list) < self.homo+1:
                # append quater seq to cannot control list
                for seq in quater_seq_list:
                    # "0123" at the start position is the mapping placeholder
                    cannot_control_list.append("0123" + seq)
                break

            # list第一个元素进行切分，
            _homo_control_list = []
            quater_seq_found = []
            quater_seq_found_idnex = []

            part_seq_list = splitSeq(quater_seq_list[0])
            for part_seq in part_seq_list:
                for i in range(1, len(quater_seq_list)):
                    # "0123" at the start position is the mapping placeholder
                    quater_seq = "0123"+quater_seq_list[i]
                    # 已找到符合条件的序列，跳过
                    if i in quater_seq_found_idnex:
                        continue

                    inserted_seq = insertSeqAndJudge(part_seq, quater_seq)
                    if inserted_seq != "":
                        # quater_seq_found.append(part_seq)
                        quater_seq_found_idnex.append(i)
                        _homo_control_list.append(inserted_seq+"0")
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

        add_num = len(self.homo_control_list[0]) - len(self.cannot_control_list[0])
        control_seq_list = []

        for quater_seq in self.cannot_control_list:
            control_seq = ""
            for i in range(add_num):
                neibhor_base = []
                if control_seq != "":
                    neibhor_base.append(control_seq[-1])

                if quater_seq[i*self.homo: (i+1)*self.homo] != "":
                    neibhor_base.append(quater_seq[i*self.homo])

                ## insert
                # add random number
                for n in range(4):
                    if str(n) not in neibhor_base:
                        control_seq += str(n)
                        control_seq += quater_seq[i*self.homo: (i+1)*self.homo]
                        break

            # add 3 at the end of control_seq, the length is different with insert seq
            control_seq = control_seq[:-1] + "3"
            control_seq_list.append(control_seq)

        self.total_seq_list = self.homo_control_list + control_seq_list

    

    def chooseMappingAndGetNtSeq(self):
        nt_seq_list = []
        for seq in self.total_seq_list:
            # chose the suitable combination
            at_cha, cg_cha = [], []
            tag = len(seq)+1
            chose_combin = []
            num_dic = {str(i): seq.count(str(i)) for i in range(4)}
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
            nt_seq = "".join([map_dic[i] for i in seq])
            nt_seq_list.append(nt_seq)

        self.nt_seq_list = nt_seq_list



    def saveNtSeq(self):
        if not self.output_path.endswith(".dna"):
            self.output_path += ".dna"

        with open(self.output_path, "w") as f:
            for seq in self.nt_seq_list:
                f.write(seq+"\n")

    def saveConfig(self):
        density = os.path.getsize(self.input_path)*8/(len(self.nt_seq_list[0])*len(self.nt_seq_list))

        config_path = self.output_path + ".config"
        with open(config_path, "w") as f:
            outline = "fileExtension\t{}\nsplitLength\t{}\nrunningLength\t{}\nindexLength\t{}\nMD5\t{}\ndensity\t{}\n".format(
                      self.input_path.split(".")[-1], self.split_len, self.homo, self.idnex_len, scan(self.input_path), density)
            f.write(outline)

    def saveConfigIdentity(self):
        pass


    def main(self):
        self.transQuaternarySeq()
        self.splitQuaternarySeq()
        self.homoControlWithSeqInsert()
        self.homoControlwithRandomInsert()
        self.chooseMappingAndGetNtSeq()
        self.saveNtSeq()
        self.saveConfig()
        self.saveConfigIdentity()



class textEncode(insertEncode):
    # data file is coded by utf-8

    def convertQuaternaryNum(self, ch_num):
        pre_key= 0
        bytes_num_dic = {192:1, 12480:2, 798912:3, 51130560:4}
        for key, value in bytes_num_dic.items():
            if ch_num >= key:
                pre_key = key
                continue

            _ch_num = ch_num - pre_key
            num_list = tool().decimal2OtherSystem(_ch_num, 4, precision=value*4)

            for i in range(value-1):
                num_list[i] = "3"
            break

        return "".join(num_list)


    def transQuaternarySeq(self):
        quater_seq = ""
        with open(self.input_path, encoding="utf-8") as f:
            for line in f:
                for ch in line:
                    ch_num = ord(ch)
                    quater_ch_num = self.convertQuaternaryNum(ch_num)
                    quater_seq += quater_ch_num

        self.quater_seq = quater_seq


class imageEncode(insertEncode):
    def transQuaternarySeq(self):
        quater_seq = ""
        precision = len(tool().decimal2OtherSystem(255, 4))

        im = Image.open(self.input_path)
        # dpi = list(im.info.get("dpi"))[0]
        pix = im.load()
        width, height = im.size
        for y in range(height):
            for x in range(width):
                r, g, b = pix[x, y]
                quater_seq += "".join(tool().decimal2OtherSystem(r, 4, precision))
                quater_seq += "".join(tool().decimal2OtherSystem(g, 4, precision))
                quater_seq += "".join(tool().decimal2OtherSystem(b, 4, precision))

        self.quater_seq = quater_seq
        self.width = width
        self.height = height

    def saveConfigIdentity(self):
        config_path = self.output_path+ ".config"

        with open(config_path, "a") as f:
            outline = "width\t{}\nheight\t{}\\n".format(self.width, self.height)
            f.write(outline)


class audioEncode(insertEncode):
    
    def transQuaternarySeq(self):
        sound = AudioSegment.from_file(self.input_path)
        samples = sound.get_array_of_samples()
        shifted_sample = list(samples)

        add_num = int(sound.max_possible_amplitude/2)
        shifted_sample = [i+add_num for i in shifted_sample]

        # unit_len
        unit_len = len(tool().decimal2OtherSystem(sound.max_possible_amplitude-1,4))
        quater_list = ["".join(tool().decimal2OtherSystem(i, 4, unit_len)) for i in shifted_sample]

        quater_seq = "".join(quater_list)

        self.quater_seq = quater_seq
        self.max = int(sound.max_possible_amplitude)
        self.sample_num = len(shifted_sample)
        self.channels = sound.channels
        self.array_type = sound.array_type
        self.frame_rate = sound.frame_rate
        self.shifted_sample = shifted_sample

    def saveConfigIdentity(self):
        config_path = self.output_path+ ".config"

        with open(config_path, "a") as f:
            outline = "maxPossibleAmplitude\t{}\nsampleNumber\t{}\nchannels\t{}\narrayType\t{}\nframeRate\t{}\n".format(
                        self.max, self.sample_num, self.channels, self.array_type, self.frame_rate)
            f.write(outline)


class videoEncode:
    def __init__(self, input_path, output_dir, split_len=200, homo=4, gc_standard = 0.5):
        file_name = os.path.basename(input_path).replace(".", "_")
        output_dir = "{}{}{}".format(output_dir, os.sep, file_name)
        tmp_dir = "{}{}{}".format(output_dir, os.sep, "temp")
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


    def getAudio(self):
        audio_path = "{}{}{}.mp3".format(self.tmp_dir, os.sep, self.file_name)
        clip = mp.VideoFileClip(self.input_path)
        audio = clip.audio
        audio.write_audiofile(audio_path)
        clip.close()

        self.audio_path = audio_path


    def getImages(self):
        image_path_list = []
        cap = cv2.VideoCapture(self.input_path)
        fps = cap.get(cv2.CAP_PROP_FPS)

        # video config
        config_path = self.output_dir + os.sep + self.file_name + ".config"
        outline = "fileExtension\t{}\nFPS\t{}\n".format(self.input_path.split(".")[-1], int(fps))
        with open(config_path, "w") as f:
            f.write(outline)


        count = 0    
        while True:
            flag, frame = cap.read()
            if flag == False:
                break

            # no Chinese
            # image_path = "{}{}{}_{}.jpg".format("D:/", os.sep, self.file_name, count)
            image_path = "{}{}{}_{}.jpg".format(self.tmp_dir, os.sep, self.file_name, count)
            image_path_list.append(image_path)
            cv2.imwrite(image_path, frame)

            count += 1

        cap.release()



        self.image_path_list = image_path_list


    def codingFiles(self):
        # audio
        audioEncode(self.audio_path, self.output_dir).main()

        # images
        for image_path in self.image_path_list:
            imageEncode(image_path, self.output_dir).main()


    def main(self):
        self.getAudio()
        self.getImages()
        self.codingFiles()



class insertDecode:
    def __init__(self, seq_dir, save_dir):
        try:
            os.makedirs(save_dir)
        except:
            pass 


        file_list = os.listdir(seq_dir)
        for file in file_list:
            if file.endswith(".dna"):
                seq_path = "{}{}{}".format(seq_dir, os.sep, file)
                save_path = "{}{}{}".format(save_dir, os.sep, file)
                break

        self.seq_path = seq_path
        self.save_path = save_path

    def readSeq(self):
        nt_seq_list = []
        with open(self.seq_path) as f:
            for seq in f:
                nt_seq_list.append(seq.strip())

        self.nt_seq_list = nt_seq_list


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



    def recoverNtSeq(self):
        remove_insert_len = 4 + int(self.config_dic["indexLength"]) + int(self.config_dic["splitLength"])
        # has 1 tag at the end of seq
        isnert_len = len(self.nt_seq_list[0])-1 - remove_insert_len

        quater_seq_list, part_seq_list = [], []
        # quater_seq_part_seq_dic = {}
        for seq in self.nt_seq_list:
            original_seq = ""
            insert_seq = ""
            tag = seq[-1]
            seq = seq[:-1]

            step = int(self.config_dic["runningLength"])+1
            for i in range(0, len(seq), step):
                unit_seq = seq[i: i+step]

                if len(insert_seq) == isnert_len:
                    original_seq += unit_seq
                    continue
                if len(original_seq) == remove_insert_len:
                    insert_seq += unit_seq
                    continue

                insert_seq += unit_seq[0]
                for nt in unit_seq[1:]:
                    if len(original_seq) == remove_insert_len:
                        insert_seq += nt 
                    else:
                        original_seq += nt 

            # trans to quaternary seq
            nt_dic = {original_seq[i]: str(i) for i in range(4)}
            quater_seq = "".join([nt_dic[i] for i in original_seq][4:])
            quater_seq_list.append(quater_seq)

            if nt_dic[tag] == "0":
                part_seq = "".join([nt_dic[i] for i in insert_seq])
                part_seq_list.append(part_seq)
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
        part_index_len = len(tool().decimal2OtherSystem(int(self.config_dic["runningLength"])-1, 4))

        for i in range(0, len(self.part_seq_list), int(self.config_dic["runningLength"])):
            index = self.part_seq_list[i][:int(self.config_dic["indexLength"])]

            unit_seq_list = self.part_seq_list[i: i+ int(self.config_dic["runningLength"])]
            index_head_unit_seq_list = [seq[-1*part_index_len:] + seq[len(index):-1*part_index_len] for seq in unit_seq_list]
            index_head_unit_seq_list.sort()

            part_recover_seq = index + "".join([seq[part_index_len:] for seq in index_head_unit_seq_list])
            part_recover_seq_list.append(part_recover_seq)

        self.part_recover_seq_list = part_recover_seq_list

    def combineSeq(self):
        total_seq_list = self.quater_seq_list + self.part_recover_seq_list
        total_seq_list.sort()
        quater_seq_complement = ""

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


    def removeTotalZeroPadding(self):
        pass



    def recoverFile(self):
        pass



    def writeFile(self):
        pass

    def md5Checkout(self):
        origin_md5 = self.config_dic["MD5"]
        decode_md5 = scan(self.save_path)

        if origin_md5 != decode_md5:
            raise ValueError("MD5 checkout failed!!\nSequence Path\t{}\nData file\t{}\nDecode file\t{}\n".format(
                self.seq_path, origin_md5, decode_md5))
        else:
            print("MD5 checkout: PASS")


    def main(self):
        self.readSeq()
        self.readConfig()
        self.recoverNtSeq()
        self.recoverPartSeq()
        self.combineSeq()
        self.removeTotalZeroPadding()
        self.recoverFile()
        self.writeFile()
        self.md5Checkout()



class textDecode(insertDecode):

    def removeTotalZeroPadding(self):
        while True:
            if self.quater_seq_complement[-1] == "0":
                self.quater_seq_complement = self.quater_seq_complement[:-1]
            else:
                break


    def recoverFile(self):
        bytes_num_dic = {0:0, 1:192, 2:12480, 3:798912, 4:51130560}

        # 连续最高位的个数决定了character的个数，4进制中最高位为3
        cha_list = []
        p1 = 0
        while True:
            if p1 == len(self.quater_seq_complement):
                break

            if int(self.quater_seq_complement[p1]) < 3:
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

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
            cha_num = self.quater_seq_complement[p2: p1+4*unit_num]
            cha_num = int(cha_num, 4) + bytes_num_dic[unit_num-1]
            cha_list.append(chr(cha_num))
            p1 += 4*unit_num

            self.cha_seq = "".join(cha_list)

    def writeFile(self):
        with open(self.save_path, "w", encoding="utf-8") as f:
            f.write(self.cha_seq)


    
class imageDecode(insertDecode):

    def recoverFile(self):
        rgb_list = []
        for i in range(0, len(self.quater_seq_complement), 4):
            num = self.quater_seq_complement[i: i+4]
            num = int(num, 4)
            rgb_list.append(num)

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

    def md5Checkout(self):
        # 恢复的图片压缩率可能不一致 但是不影响信息传递  此处不进行校验
        pass


class audioDecode(insertDecode):
    def recoverFile(self):
        max_possible_amplitude = int(self.config_dic["maxPossibleAmplitude"])
        sample_num = int(self.config_dic["sampleNumber"])

        # from zero, -1
        unit_len = len(tool().decimal2OtherSystem(max_possible_amplitude-1, 4))
        
        sample_list = []
        for i in range(sample_num):
            unit = self.quater_seq_complement[i*unit_len: (i+1)*unit_len]
            sample_list.append(int(unit, 4))

        # to array
        dtype = "int{}".format(int(math.log(max_possible_amplitude,2)+1))
        sample_array = np.array(sample_list, dtype = dtype)
        sample_array -= int(max_possible_amplitude/2)
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

        self.sound = sound


    def md5Checkout(self):
        pass


class videoDecode:

    def __init__(self, seq_dir, save_dir):
        tmp_dir = "{}{}{}_decodeTmp".format(seq_dir, os.sep, os.path.basename(seq_dir))
        try:
            os.makedirs(tmp_dir)
        except:
            pass 


        self.seq_dir = seq_dir
        self.save_path = save_dir + os.sep + os.path.basename(seq_dir) + ".mp4"
        self.tmp_dir = tmp_dir

    def recoverAudio(self):
        for file in os.listdir(self.seq_dir):
            if file.endswith("_mp3"):
                audio_dir = "{}{}{}".format(self.seq_dir, os.sep, file)
                audioDecode(audio_dir, self.tmp_dir).main()

        self.audio_path = ["{}/{}".format(self.tmp_dir, i) for i in os.listdir(self.tmp_dir) if i.endswith("mp3")][0]

    def recoverImages(self):
        for file in os.listdir(self.seq_dir):
            if file.endswith("_jpg"):
                image_dir = "{}{}{}".format(self.seq_dir, os.sep, file)
                imageDecode(image_dir, self.tmp_dir).main()

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
        video = mp.VideoFileClip(self.save_path)
        audio = mp.AudioFileClip(self.audio_path)
        video = video.set_audio(audio)

        video.write_videofile(self.save_path)




    def main(self):
        self.recoverAudio()
        self.recoverImages()
        self.recoverVideo()
        self.combineVideoAudio()


def encode(input_path, encode_dir, split_len=200, homo=4, gc_standard = 0.5):
    file_extension = os.path.basename(input_path).split(".")[-1]

    if file_extension in ["txt"]:
        cl = textEncode(input_path, encode_dir, split_len = split_len, homo = homo, gc_standard = gc_standard)
    elif file_extension in ["jpg", "jpeg", "png"]:
        cl = imageEncode(input_path, encode_dir, split_len = split_len, homo = homo, gc_standard = gc_standard)
    elif file_extension in ["mp3", "wav"]:
        cl = audioEncode(input_path, encode_dir, split_len = split_len, homo = homo, gc_standard = gc_standard)
    elif file_extension in ["mp4"]:
        cl = videoEncode(input_path, encode_dir, split_len = split_len, homo = homo, gc_standard = gc_standard)
    else:
        raise ValueError("The type of {} file can't be processed!".format(file_extension))

    cl.main()


def decode(encode_dir, decode_dir):
    _cl = insertDecode(encode_dir, decode_dir)
    _cl.readConfig()
    file_extension = _cl.config_dic["fileExtension"]

    if file_extension in ["txt"]:
        cl = textEncode(encode_dir, decode_dir)
    elif file_extension in ["jpg", "jpeg", "png"]:
        cl = imageEncode(encode_dir, decode_dir)
    elif file_extension in ["mp3", "wav"]:
        cl = audioEncode(encode_dir, decode_dir)
    elif file_extension in ["mp4"]:
        cl = videoEncode(encode_dir, decode_dir)
    else:
        raise ValueError("The config file is not correct!\nFile extension: {}".format(file_extension))

    cl.main()




if __name__ == "__main__":
    input_path = "D:/BaiduSyncdisk/中科院先进院合成所/项目/2021_12_06-分割插入编码/test_data/English.txt"
    encode_dir = os.path.dirname(input_path)
    decode_dir = encode_dir+os.sep+os.path.basename(input_path).replace(".", "_")
    # # encode_path = ".".join(input_path.split(".")[:-1])+".dna"
    # # decode_path = encode_path + "."+input_path.split(".")[-1]

    cl = textEncode(input_path, encode_dir)
    cl.main()

    # quater_seq=cl.quater_seq
    # quater_seq_list = cl.quater_seq_list
    # homo_control_list = cl.homo_control_list
    # cannot_control_list = cl.cannot_control_list
    # total_seq_list = cl.total_seq_list
    # nt_seq_list = cl.nt_seq_list
    # shifted_sample = cl.shifted_sample

    # gc_list = [(i.count("C")+i.count("G"))/len(i) for i in nt_seq_list]
    # gc_list.sort()
    # print(round(gc_list[0], 4), round(gc_list[-1],4))

    #########################  decode
    encode_dir = os.path.dirname(input_path)+os.sep+os.path.basename(input_path).replace(".", "_")
    decode_dir = encode_dir

    cld = textDecode(encode_dir, decode_dir)
    cld.main()

    # quater_seq_list_d = cld.quater_seq_list 
    # part_seq_list = cld.part_seq_list
    # part_recover_seq_list = cld.part_recover_seq_list
    # total_seq_list_d = cld.total_seq_list
    # quater_seq_complement = cld.quater_seq_complement
    # nt_seq_list_d = cld.nt_seq_list
    # sound = cld.sound



