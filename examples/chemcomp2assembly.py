#!/usr/bin/python -u
#
#

import os
import sys
import json

import pprint

#
# input is a list (top-level json object) ot lists: saveframe categories.
# there should be only one sublist w/ only one dict in it: the single NMR-STAR
#   chem. comp. created from a MOL/SDF.
#
# this adds a single entity saveframe and an assembly saveframe with a few 
#   boilerplate values. Atom list is copied from chem. comp. to assembly
#   b/c assembly_atom_id is technically the parent key for all atoms in the
#   entry so we populate it.
#

if __name__ == "__main__" :

    dat = json.load( sys.stdin )

    comp = dat[0][0]

# fix up if not set:
#
    if not "aromatic" in comp.keys() :
        arom = False 
        for a in comp["atoms"] :
            if a["aromatic"] :
                arom = True
                break
# JIC
#
        if not arom :
            for b in comp["bonds"] :
                if b["arom"] == "yes" :
                    arom = True
                    break

        comp["aromatic"] = arom

# entity: boilerplate
#
    enty = {}
    enty["formula_weight"] = comp["weight"]
    enty["nonpolymer_comp_label"] = comp["name"]
    enty["paramagnetic"] = comp["paramagnetic"]
    enty["id"] = 1
    enty["name"] = "entity_1"
    enty["sf_category"] = "entity" 
    enty["thiol_state"] = "not present"
    enty["type"] = "non-polymer"
    enty["ambiguous_conformational_states"] = "no"
    enty["ambiguous_chem_comp_sites"] = "no"
    enty["nstd_linkage"]  = "no"
    enty["number_of_monomers"] = 1
    enty["number_of_nonpolymer_components"] = 1
    enty["entity_comp_index"] = [ {
        "auth_seq_id": 1,
        "comp_label": comp["name"],
        "entity_id": 1,
        "id": 1
    } ]

    dat.insert( 0, [enty] )

# assembly: copy of atoms
#

    assy = {}
    assy["id"] = 1
    assy["name"] = "assembly_1"
    assy["sf_category"] = "assembly"
    assy["number_of_components"] = 1
    assy["paramagnetic"] = enty["paramagnetic"]
    assy["thiol_state"] = enty["thiol_state"]
    assy["physical_state"] = "native"
    assy["chemical_exchange_state"] = "no"
    assy["conformational_isomer"] = "no"
    assy["atoms"] = []
    for a in comp["atoms"] :
        assy["atoms"].append( {
            "assembly_atom_id": a["ordinal"],
            "id": a["id"], 
            "type_symbol": a["type"],
            "assembly_id": 1, 
            "comp_index_id": 1, 
            "entity_assembly_id": 1, 
            "entity_id": 1, 
            "seq_id": 1, 
        } )

    dat.insert( 0, [enty] )

    json.dump( dat, sys.stdout, indent = 2, sort_keys = True )


#
#
#
