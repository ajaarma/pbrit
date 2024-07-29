#!/usr/bin/python

import re
###############################################
#
# Class ONTOLOGY
#
#
# Description: Generic ONTOLOGY class providing methods to process specific data structure of OBO files which includes GO, HPO, DO and MPO
#
#
#
# Things to do: (1) Improvemnet towards more OOP style programming.
#               (2) Try to create single method that  can handle processing of GO, HPO, MPO, DO all together. Currently there are several repetitions of the code.
# 
##############################################

class OBOEngine:

    def __init__(self, elements=[]):

        self.__elements={}
        for e in elements:
            self.__elements[e]=1

    def getGOHashTable(self,obo_file):

        self.obo_file = obo_file
        self.GO_hash = {}

        go_id = 0
        namespaces=names=defs=is_a=''
        is_a_str=[]=consider_obj=alt_id_obj=rel_part=is_obs_str=rep_by_str=rel_regulate=[]
        id_pattern = re.compile('^id\:\W+GO\:(\w+)$')
        namespace_pattern = re.compile('^namespace\:\W+(.+)$')
        name_pattern = re.compile('^name\:\W+(.+)$')
        def_pattern = re.compile('^def\:\W+(.+)$')
        is_a_pattern = re.compile('^is_a\:\W+GO\:(\w+)\W+\!\W+(\w+)')
        is_obs_pattern = re.compile('^is_obsolete\:\W+true')
        is_rep_pattern = re.compile('^replaced_by\:\W+GO\:(\w+)\W*$')
        is_consider_pattern = re.compile('^consider\:\W+GO\:(\w+)\W*$')
        is_alt_pattern = re.compile('^alt_id\:\W+GO\:(\w+)\W*')
        is_relp_pattern = re.compile('^relationship\:\W+part\_of\W+GO\:(\w+)\W*\!')
        is_relg_pattern = re.compile('^relationship\:\W+regulates\W+GO\:(\w+)\W*\!')

        try:
            fh = open(self.obo_file)

            for lines in fh:
                lines = lines.strip()
                exp = re.search('\W+',lines)
                if exp:
                    if id_pattern.search(lines):
                        strs = re.split('\:(\w+)$',lines)
                        go_id = strs[1].strip()
                        continue

                    if namespace_pattern.search(lines):
                        strs = re.split('\:\W+(\w+)$',lines)
                        namespaces = strs[1].strip()
                        continue

                    if name_pattern.search(lines):  
                        strs = re.split('\:\W+(.+)$',lines)
                        names=strs[1].strip()
                        continue

                    if def_pattern.search(lines):
                        strs = re.split('\:\W+(.+)$',lines)
                        defs = strs[1].strip()
                        continue
                    if is_a_pattern.search(lines):
                        strs = re.split('^is_a\:\W+GO\:(\w+)\W+(.+)',lines,2)
                        is_a = strs[1].strip()
                        is_a_str.append(is_a)
                        continue
                    if is_obs_pattern.search(lines):
                        strs = re.split('\:\W+',lines)
                        is_obs = strs[1]
                        is_obs_str.append(is_obs)
                        continue
                    if is_rep_pattern.search(lines):
                        strs = re.split('^replaced_by\:\W+GO\:(\w+)\W*',lines)
                        replaced_by = strs[1].strip()
                        rep_by_str.append(replaced_by)
                        continue
                    if is_consider_pattern.search(lines):
                        strs = re.split('^consider\:\W+GO\:(\w+)\W*',lines)
                        consider=strs[1].strip()
                        consider_obj.append(consider)
                        continue
                    if is_alt_pattern.search(lines):
                        strs = re.split('^alt_id\:\W+GO\:(\w+)\W*',lines)
                        alt_id = strs[1].strip()
                        alt_id_obj.append(alt_id)
                        continue
                    if is_relp_pattern.search(lines):
                        strs = re.split('^relationship\:\W+part\_of\W+GO\:(\w+)\W*\!',lines)
                        relationship=strs[1].strip()
                        rel_part.append(relationship)
                        continue
                    if is_relg_pattern.search(lines):
                        strs =re.split('^relationship\:\W+regulates\W+GO\:(\w+)\W*\!',lines)
                        regulate=strs[1].strip()
                        rel_regulate.append(regulate)
                        continue
                    else:
                        if go_id:
                            self.GO_hash[go_id]={'namespaces':namespaces,
                                                 'names':names,'def':defs,
                                                 'is_a':is_a_str,
                                                 'consider':consider_obj,
                                                 'alt_id':alt_id_obj,
                                                 'part_of':rel_part,
                                                 'regulate':rel_regulate,
                                                 'is_obsolete':is_obs_str,
                                                 'replaced_by':rep_by_str}
                            is_a_str = []
                            consider_obj = []
                            alt_id_obj=[]
                            rel_part=[]
                            rel_regulate=[]
                            is_obs_str=[]
                            rep_by_str=[]
                            continue

        except IOError:
            import sys
            print "No such file exists: ",obo_file
            print "Please enter a valid OBO file "
        else:
            fh.close()

        return self.GO_hash



    def getHPOHashTable(self,obo_file):
        self.obo_file = obo_file
        self.HP_hash = {}

        hp_id = 0
        names=defs=is_a=''
        is_a_str=consider_obj=alt_id_obj=rel_part=is_obs_str=rep_by_str=rel_regulate=[]

        id_pattern = re.compile('^id\:\W+HP\:(\w+)$')
        name_pattern  = re.compile('^name\:\W+(.+)$')
        def_pattern = re.compile('^def\:\W+(.+)$')
        is_a_pattern = re.compile('^is_a\:\W+HP\:(\w+)\W+\!\W+(\w+)')
        is_obs_pattern = re.compile('^is_obsolete\:\W+true')
        is_alt_pattern = re.compile('^alt_id\:\W+HP\:(\w+)\W*')

        try:
            fh = open(self.obo_file)
            for lines in fh:
                lines = lines.strip()
                exp = re.search('\W+',lines)

                if exp:
                    if id_pattern.search(lines):
                        strs = re.split('\:(\w+)$',lines)
                        hp_id = strs[1].strip()
                        continue
                    if name_pattern.search(lines):
                        strs = re.split('\:\W+(.+)$',lines)
                        names = strs[1].strip()
                        continue
                    if def_pattern.search(lines):
                        strs = re.split('\:\W+(.+)$',lines)
                        defs = strs[1].strip()
                        continue
                    if is_a_pattern.search(lines):
                        strs = re.split('^is_a\:\W+HP\:(\w+)\W+(.+)',lines,2)
                        is_a = strs[1]
                        is_a_str.append(is_a)
                        continue
                    if is_obs_pattern.search(lines):
                        strs = re.split('\:\W+',lines)
                        is_obs = strs[1].strip()
                        is_obs_str.append(is_obs)
                        continue
                    if is_alt_pattern.search(lines):
                        strs = re.split('^alt_id\:\W+HP\:(\w+)\W*',lines)
                        alt_id = strs[1].strip()
                        alt_id_obj.append(alt_id)
                        continue
                else:
                    if hp_id:
                        self.HP_hash[hp_id]={'names':names,
                                           'def':defs, 
                                           'is_a':is_a_str,
                                           'alt_id':alt_id_obj,
                                           'is_obsolete':is_obs_str}

                        is_a_str = []
                        alt_id_obj = []
                        is_obs_str = []
                        names = []
                        defs = []
                        continue

        except (IOError, TypeError):
            import sys
            print "NO such file exists: ",self.obo_file
            print "Please enter a valid OBO file"
            sys.exit()
        else:
            fh.close()

        return self.HP_hash

    def getDOHashTable(self,obo_file):

        self.obo_file = obo_file
        self.DO_hash = {}

        do_id = 0
        names=defs=is_a=''
        is_a_str=consider_obj=alt_id_obj=rel_part=is_obs_str=rep_by_str=rel_regulate=[]

        id_pattern = re.compile('^id\:\W+DOID\:(\w+)$')
        name_pattern  = re.compile('^name\:\W+(.+)$')
        def_pattern = re.compile('^def\:\W+(.+)$')
        is_a_pattern = re.compile('^is_a\:\W+DOID\:(\w+)\W+\!\W+(\w+)')
        is_obs_pattern = re.compile('^is_obsolete\:\W+true')
        is_alt_pattern = re.compile('^alt_id\:\W+DOID\:(\w+)\W*')

        try:
            fh = open(self.obo_file)
            for lines in fh:
                lines = lines.strip()
                exp = re.search('\W+',lines)

                if exp:
                    if id_pattern.search(lines):
                        strs = re.split('\:(\w+)$',lines)
                        do_id = strs[1].strip()
                        continue
                    if name_pattern.search(lines):
                        strs = re.split('\:\W+(.+)$',lines)
                        names = strs[1].strip()
                        continue
                    if def_pattern.search(lines):
                        strs = re.split('\:\W+(.+)$',lines)
                        defs = strs[1].strip()
                        continue
                    if is_a_pattern.search(lines):
                        strs = re.split('^is_a\:\W+DOID\:(\w+)\W+(.+)',lines,2)
                        is_a = strs[1]
                        is_a_str.append(is_a)
                        continue
                    if is_obs_pattern.search(lines):
                        strs = re.split('\:\W+',lines)
                        is_obs = strs[1].strip()
                        is_obs_str.append(is_obs)
                        continue
                    if is_alt_pattern.search(lines):
                        strs = re.split('^alt_id\:\W+DOID\:(\w+)\W*',lines)
                        alt_id = strs[1].strip()
                        alt_id_obj.append(alt_id)
                        continue
                else:
                    if do_id:
                        self.DO_hash[do_id]={'names':names,
                                           'def':defs, 
                                           'is_a':is_a_str,
                                           'alt_id':alt_id_obj,
                                           'is_obsolete':is_obs_str}

                        is_a_str = []
                        alt_id_obj = []
                        is_obs_str = []
                        names = []
                        defs = []
                        continue



        except (IOError, TypeError):
            import sys
            print "NO such file exists: ",self.obo_file
            print "Please enter a valid OBO file"
            sys.exit()
        else:
            fh.close()

        return self.DO_hash


    def getMPOHashTable(self,obo_file):

        self.obo_file = obo_file
        self.MP_hash = {}

        mp_id = 0
        names=defs=is_a=''
        is_a_str=consider_obj=alt_id_obj=rel_part=is_obs_str=rep_by_str=rel_regulate=[]

        id_pattern = re.compile('^id\:\W+MP\:(\w+)$')
        name_pattern  = re.compile('^name\:\W+(.+)$')
        def_pattern = re.compile('^def\:\W+(.+)$')
        is_a_pattern = re.compile('^is_a\:\W+MP\:(\w+)\W+\!\W+(\w+)')
        is_obs_pattern = re.compile('^is_obsolete\:\W+true')
        is_alt_pattern = re.compile('^alt_id\:\W+MP\:(\w+)\W*')

        try:
            fh = open(self.obo_file)
            for lines in fh:
                lines = lines.strip()
                exp = re.search('\W+',lines)

                if exp:
                    if id_pattern.search(lines):
                        strs = re.split('\:(\w+)$',lines)
                        mp_id = strs[1].strip()
                        continue
                    if name_pattern.search(lines):
                        strs = re.split('\:\W+(.+)$',lines)
                        names = strs[1].strip()
                        continue
                    if def_pattern.search(lines):
                        strs = re.split('\:\W+(.+)$',lines)
                        defs = strs[1].strip()
                        continue
                    if is_a_pattern.search(lines):
                        strs = re.split('^is_a\:\W+MP\:(\w+)\W+(.+)',lines,2)
                        is_a = strs[1]
                        is_a_str.append(is_a)
                        continue
                    if is_obs_pattern.search(lines):
                        strs = re.split('\:\W+',lines)
                        is_obs = strs[1].strip()
                        is_obs_str.append(is_obs)
                        continue
                    if is_alt_pattern.search(lines):
                        strs = re.split('^alt_id\:\W+MP\:(\w+)\W*',lines)
                        alt_id = strs[1].strip()
                        alt_id_obj.append(alt_id)
                        continue
                else:
                    if mp_id:
                        self.MP_hash[mp_id]={'names':names,
                                           'def':defs, 
                                           'is_a':is_a_str,
                                           'alt_id':alt_id_obj,
                                           'is_obsolete':is_obs_str}

                        is_a_str = []
                        alt_id_obj = []
                        is_obs_str = []
                        names = []
                        defs = []
                        continue

        except (IOError, TypeError):
            import sys
            print "NO such file exists: ",self.obo_file
            print "Please enter a valid OBO file"
            sys.exit()
        else:
            fh.close()

        return self.MP_hash

