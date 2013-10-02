#! /usr/bin/env python

# File: newick.py
# Author: Dylan Schwilk (www.pricklysoft.org)
# Date: 2005/03/09
# Copyright 2003, 2004, 2005 Dylan W. Schwilk

# GNU
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version.
#   
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

"""Module for reading newick - formatted trees.

   Functions:
   
   create_tree(src) - will accept a token
                     list, a string or a file-like object to read
                     from.
                     
   read_trees(src) - reads list of trees separated by semi-colons as in a *new
                     file.
"""

__author__  =    '''Dylan Schwilk (www.pricklysoft.org)'''

from shlex import shlex
from types import StringType, ListType
from cStringIO import StringIO
from phylotree import PhyloTree


def read_trees(src) :
    trees = []
    strs = src.split(";")
    for s in strs:
        s = s.strip()
        if len(s) > 2:  # skip empty lines
            trees.append(create_tree(s))
    return trees

def create_tree(l) :
    '''Reads Newick format tree from token list, string, or file-like object
       The function does not check for comments and expects an already-cleaned
       up tree.'''

    if not type(l) is ListType :
        l = parse(l)
    
    root = PhyloTree()
    node = root
    lp = 0
    rp = 0
    t = 0
    while t < len(l) :
        if l[t] == '(' :
            lp += 1
            newnode = PhyloTree()
            node.add_child(newnode)
            node = newnode
        elif l[t] == ')' :
            rp += 1
            node = node.parent
        elif l[t] == ',' :
            newnode = PhyloTree()
            node.parent.add_child(newnode)
            node = newnode   
        elif l[t] == ':' :
            t += 1
            node.bl =  float(l[t])
        else :
            node.label = l[t]
            #print node.label
        t += 1
    if lp != rp :
        raise StandardError('Unbalanced parentheses in tree description',
                            ''.join(map(str,l)))
    return root
                

# ------- Private functions ------------ #
class Tokenizer(shlex):
    """Provides tokens for parsing Newick-format trees"""
    def __init__(self, infile):
        shlex.__init__(self, infile)
        self.commenters = ''
        self.wordchars = self.wordchars+'-.'
        self.quotes = "'"

def parse(src):
    """
    Parse a Newick-formatted tree description
    input is any file-like object that can be coerced into shlex,
    or a string (converted to StringIO)
    """
    if type(input) is StringType:
        src = StringIO(src)
    
    # start_pos = src.tell()
    tokens = Tokenizer(src)
    result = []
    while 1:
        token = tokens.get_token()
        if token == '' : break
        result.append(token)
    return result
        

# Main Test function               
if __name__ == '__main__':
    from phylotree import equivalent

    print "Phylotree test\n"

    treestr = "(Ephedra_triandra,(((Zephyranthes_longystila,Nothoscordum_inodorum,Hypoxis_humilis),((Sisyrinchium_chilensis,Sisyrinchium_unguiculatum),(((Deyeuxia_hieronymi,(Festuca_hieronymi,Festuca_tucumanica),Cortaderia_rudiuscula,Agrostis_montevidensis,(Piptochaetium_montevidense,Piptochaetium_chico),Sorghum_halepense,(Setaria_geniculata,Setaria_pampeana),Bothriochloa_laguroides,Briza_subaristata,Neobouteloua_lophostachya,Tripogon_spicatum,(Paspalum_notatum,Paspalum_quadrifalium),Vulpia_myurus,Monanthochloe_litoralis,Pappophorum_chaco,(Aristida_achalensis,Aristida_adscencionis,Aristida_mendocina,Aristida_spegazzini),Chloris_retusa,Muhlenbergia_peruviana,Axonopus_indet,(Schizachyrium_microstachyum,Schizachyrium_gris,Schizachyrium_verged,Schizachyrium_tenerum),Melica_macra,Gouinia_paraguayensis,(Bromus_auleticus,Bromus_brevis,Bromus_catharticus),Sorghastrum_pellitum,(Stipa_eriostachya,Stipa_filiculmis,Stipa_flexibarbata,Stipa_neesiana,Stipa_tenuissima,Stipa_trichotoma),Bouteloa_aristidoides,Cenchrus_pauciflorus,(Eragrostis_lugens,Eragrostis_peludo),(Poa_resinulosa,Poa_stuckertii),(Sporobolus_indicus,Sporobolus_piramidatus),Trichloris_crinita),((Eleocharis_albibracteata,Bulbostylis_juncoides,Carex_fuscula,Cyperus_reflexus),(Juncus_achalensis,Juncus_uruguensis)),(Deutheroconia_longipetala,(Tillandsia_brioides,Tillandsia_duratii,Tillandsia_nueva),Deynacanthon_urbanianum)),Commelina_virginica))),((Rumex_acetosella,(((Suaeda_divaricata,Iresine_difusa,Heterostachys_ritteriana,(Alternanthera_pumilla,Alternanthera_pungens),Gomphrena_pulchella,Pfaffia_gnaphaloides,Chenopodium_album,Allenrolfea_patagonica,Atriplex_argentina),(Cardionema_ramosissimun,Spergula_ramosa,Cerastium_arvense)),(Allionia_incarnata,(((Portulaca_grandiflora,Portulaca_oleracea),Grahamia_bracteata,Talinum_polygaloides),(Stetsonia_coryne,Cereus_validus,Opuntia_sulfurea,Gymnocalycium_rosa,Tephrocactus_glomeratus))))),Ligaria_cuneifolia,((Oenothera_indecora,Cuphea_glutinosa),((Bulnesia_retama,Tribulus_terrestris,(Larrea_cuneifolia,Larrea_divaricata),Plectocarpa_tetracantha),((Maytenus_spinosa,Maytenus_vitisideae),Oxallis_sexenata,(Cordobia_argentea,(Euphorbia_serpens,((Croton_argentinum,Croton_sarcopetalus),Acalypha_communis)))),(((Trifolium_amabile,Trifolium_repens),Stylosanthes_gracilis,(Acacia_aroma,Acacia_caven,Acacia_fucatispina),Mimozyganthus_carinatus,(Prosopis_algarrobilla,Prosopis_flexuosa,Prosopis_torquata),Cercidium_australe,Zuccagnia_punctata,Vicia_graminea,Astragalus_parodii,Rynchosia_senna,Cologamia_ovalifolia,Senna_aphylla,Geoffroea_decorticans,(Adesmia_corimbosa,Adesmia_repetidora),Galactia_marginalis),((Margyricarpus_pinnatus,Alchemilla_pinnata),((Colletia_spinosissima,Condalia_microphylla,Ziziphus_mistol),(Celtis_pallida,Celtis_tala))))),((Lepidium_bonariensis,Capparis_atamisquea),(Tricomaria_usillo,Pavonia_aurigloba,(Sida_indet,Sida_prostrata)),(Fagara_coco,(Lithraea_terniifolia,Schinopsis_haenkeana)))),(((((Mitracarpus_cuspidatus,Mitracarpus_indet),(Borreria_verticillata,Borreria_eryngioides),Relbunium_richardianum),(Gentiana_parviflora,Aspidosperma_quebrachoblanco)),((Pithecoctenium_synanchoides,Lepechinia_floribunda,((Glandularia_dissecta,Glandularia_peruviana),Lippia_salsa,Aloysia_gratisima),(Stenandrium_dulce,Justicia_squarrosa),(Plantago_australis,Plantago_brasiliensis)),((Lycium_elongatum,Solanum_eleagnifolium,Nierenbergia_hippomanica),(Dichondra_serisea,Ipomoea_azul,Evolvulus_sericeous)))),(((Eryngium_agavifolium,Eryngium_horridum,Eryngium_nudicaule),Oreomyrrhis_andicola),(Acicarpha_tribuloides,(Carduus_thoermeri,Stevia_satureiifolia,Flourensia_campestris,Gnaphalium_gaudichaudianum,Schkuhria_pinnata,Hieracium_giganteum,Chevreulia_sarmentosa,Taraxacum_officinale,Flaveria_bidentis,(Baccharis_articulata,Baccharis_coridifolia,Baccharis_rufescens),Solidago_chilensis,(Noticastrum_argenteum,Noticastrum_marginatum),Zinnia_peruviana,Cotula_mexicana,Conyza_indet,Vernonia_nudiflora,Parthenium_hysterophorum,Bidens_triplinervia,Cyclolepis_genistoides,Heterothalamus_allienus,Spilanthes_decumbens,Chaptalia_integerrima,(Tagetes_minuta,Tagetes_argentina),Verbesina_encelioides,(Hypochoeris_argentina,Hypochoeris_caespitosa),Eupatorium_buniifolium,Achyrocline_indet,(Gamochaeta_discolor,Gamochaeta_nodiscolor),Ambrosia_tenuifolia,Lucilia_acutifolia,Ophryosporus_axilliflorus)))))))"

    print "Tree read/write test ...",
    tree = create_tree(treestr)
    if create_tree(tree.write(True)).write(False) == treestr : print "Passed"
    else : print "fail"
 

    print "Prune test ... ",
    prune_list = ["Sporobolus_indicus","Sporobolus_piramidatus","Paspalum_notatum", "Lithraea_terniifolia"]
    tree.prune_taxa(prune_list)
    if len(tree.leaves()) == 207 and len(tree.descendants())== 296: print "passed"
    else : print "fail"
    
    print 'copy, reverse and test for equivalence ....',
    tree2 = tree.copy()
    tree2.reverse()
    if equivalent(tree,tree2, True) : print "passed"
    else : print "fail"


  #  print "Reroot test ...",
  #  tree3 = tree.descendants()[18].reroot()
  #  print len(tree3.leaves()), len(tree3.descendants())
