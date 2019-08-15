#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

from __future__ import print_function
import argparse
import sys
import os
import re

parser = argparse.ArgumentParser(
    usage='%(prog)s --file_type [gff|gtf|saf] --in_file <path/to/your/db_file>',
    description="This script creates three column file: "
                "Gene ID, Gene Name and Gene Biotype",
    add_help=True
)
parser.add_argument(
    '-i',
    '--in_file',
    required=True,
    metavar='FILE',
    help="path to your annotation file"
)
parser.add_argument(
    '-t',
    '--file_type',
    metavar='STR',
    required=True,
    help="specify file type [gff|gtf|saf]"
)

parser.add_argument(
    '-f',
    '--feat_types',
    type=str,
    dest='feat_types',
    metavar='FEATURE_NAME',
    nargs='+',
    default=['gene', 'exon', 'transcript', 'CDS', 'mRNA'],
    help="specify feature types (space seperated, eg -f gene exon CDS transcript). This is third column in the input annotation file. Default [gene exon transcript CDS]"
)

def get_gtf(handler, feat_types=None, name_attr='gene_name'):

    if feat_types is None:
        feat_types = ['gene', 'exon', 'transcript', 'CDS']

    genes_attr = {}

    typesRegex = "\s.([A-z0-9 _., \- /\(\)]+)"

    for key, value in list(biotypes.items()):
        tweak = value + typesRegex
        biotypes[key] = tweak

    for i in handler:
        line = i.strip()

        if line.startswith("#"):
            continue

        feature = line.split('\t')
        if len(feature) == 9:
            if feature[2] in feat_types:
                ninthField = feature[8]

                chk_gene_id = re.search('gene_id\s.([A-z0-9_ \. \- \(\)]+)',
                                        ninthField)
                chk_gene_name = re.search('%s\s.([A-z0-9_ \. \: \- \(\)]+)' % name_attr,
                                          ninthField)

                gene_name = 'NA'
                gene_biotype = 'NA'
                chrom = feature[0]

                if chk_gene_name:
                    gene_name = chk_gene_name.group(1)

                for value in list(biotypes.values()):
                    checkBiotype = re.search(value, ninthField)
                    if checkBiotype:
                        gene_biotype = checkBiotype.group(1)

                if chrom not in genes_attr:
                    genes_attr[chrom] = {}

                if chk_gene_id:
                    gene_id = chk_gene_id.group(1)

                    if gene_id not in genes_attr[chrom]:
                        genes_attr[chrom][gene_id] = {}
                    # TODO I feel that this is error prone..
                    # I'mm imaginig a situations where previous gene_name had
                    # proper value but in the next iteration (next line) gene_name
                    # tag doesn't exists and therefor gets NA value and thats now
                    # overides proper values already stored in the dictionary
                    # in theory every line in gtf files is self contained and this
                    # shouldn't happened but my custome gtf violated that rule,
                    # but former lines would always hold "proper" value, so I guess
                    # I'll leave it as is for now
                    genes_attr[chrom][gene_id]["gene_name"] = gene_name
                    genes_attr[chrom][gene_id]["biotype"] = gene_biotype
        else:
            print("WARNING: Length of the line %d, should be 9" % len(feature), file = sys.stderr)

    return genes_attr


def create_gff_tree(handler):
    """
    Reads a GFF file handle and returns a dictionary like:
    
      {id, {biotype, parent}}
    
    The biotype or parent may be None.

    'id' is the ID attribute from the GFF, 'biotype' is the gene_biotype attribute,
    or in the case where there is no gene_biotype, the attribute is taken from the Parent (or Parent's Parent),
    walking up the hierarchy of nested features (only two branches are traversed).

    For an NCBI RefSeq type GFF, this effectively propagates biotype from gene -> mRNA -> CDS,
    (even though the CDS does not have a biotype attribute itself).
    """
    biotype_parent_tree = {}

    for i in handler:
        line = i.strip()

        if line.startswith('#'):
            continue

        fields = line.split('\t')
        contig = fields[0]
        feature_type = fields[2]
        attributes = [tuple(a.split('=')) for a in fields[8].split(';')]
        attributes = dict(attributes)

        # print(feature_type, attributes)

        ID = attributes.get('ID', None)
        biotype = attributes.get('gene_biotype', None)
        parent = attributes.get('Parent', None)
        if ID and ID not in biotype_parent_tree:
            biotype_parent_tree[ID] = {'parent': parent, 'biotype': biotype}     

    # HACK: This would be smarter as a recursive function (with some depth limit).
    #       Instead we just run this twice to propagate biotype from gene -> mRNA -> CDS
    #       for NCBI RefSeq type GFFs.
    for ID, link in list(biotype_parent_tree.items()):
        parent = link['parent']
        if not link['biotype'] and parent:
            biotype_parent_tree[ID]['biotype'] = biotype_parent_tree[parent]['biotype']

    for ID, link in list(biotype_parent_tree.items()):
        parent = link['parent']
        if not link['biotype'] and parent:
            biotype_parent_tree[ID]['biotype'] = biotype_parent_tree[parent]['biotype']

    return biotype_parent_tree


def get_gff(handler, feat_types=None, biotype_id_mapping=None):

    if feat_types is None:
        feat_types = ['gene', 'exon', 'transcript', 'CDS', 'mRNA']

    genes_attr = {}

    typesRegex = "=([A-z0-9 _., \- /\(\)]+)"

    for key, value in list(biotypes.items()):
        tweak = value + typesRegex
        biotypes[key] = tweak

    for i in handler:
        line = i.strip()

        if line.startswith('#'):
            continue

        feature = line.split('\t')
        if len(feature) == 9:
            if feature[2] in feat_types:
                ninthField = feature[8]
    
                chk_gene_id = re.search('ID=([A-z0-9_ \. \- \(\)]+)', ninthField)
                chk_gene_name = re.search('Name=([A-z0-9_ \. \: \- \(\)]+)', ninthField)
    
                gene_name = 'NA'
                gene_biotype = 'NA'
                chrom = feature[0]

                if chk_gene_name:
                    gene_name = chk_gene_name.group(1)
                
                if chrom not in genes_attr:
                    genes_attr[chrom] = {}

                if chk_gene_id:
                    gene_id = chk_gene_id.group(1)

                    if gene_id not in genes_attr[chrom]:
                        genes_attr[chrom][gene_id] = {}

                    genes_attr[chrom][gene_id]["gene_name"] = gene_name
                
                if biotype_id_mapping:
                    gene_biotype = biotype_id_mapping.get(gene_id, {}).get('biotype', 'NA')
                else:
                    for value in list(biotypes.values()):
                        checkBiotype = re.search(value, ninthField)
                        if checkBiotype:
                            gene_biotype = checkBiotype.group(1)
                
                genes_attr[chrom][gene_id]["biotype"] = gene_biotype

                # if gene_id:
                #    if gene_id.group(1) not in genes_attr:
                #        genes_attr[gene_id.group(1)] = []
                #        genes_attr[gene_id.group(1)].append(gene_name)
                #        genes_attr[gene_id.group(1)].append(gene_biotype)
                #        genes_attr[gene_id.group(1)].append(feature[0])
        else:
            print("WARNING: Length of the line %d, should be 9" % len(feature), file = sys.stderr)

    return genes_attr


def get_saf(handler):
    genes_attr = {}

    for i in handler:
        line = i.strip()
        if line.startswith('#'):
            continue

        line = line.split("\t")
        gene_id = line[0]
        chrom = line[1]
        gene_name = line[5]
        biotype = line[6]

        genes_attr[chrom][gene_id] = {"gene_name": gene_name,
                                      "biotype": biotype}

    return genes_attr


args = parser.parse_args()
file_type = args.file_type
in_file = args.in_file
feat_types = args.feat_types

biotypes = {
    "gb": 'gene_biotype',
    "gt": 'gene_type',
    "tb": 'transcript_biotype',
    "tt": 'transcript_type'
}

biotype_id_mapping = None
if file_type == "gff":
    with open(in_file) as fh:
        biotype_id_mapping = create_gff_tree(fh)

    #for k, v in biotype_id_mapping.items():
    #    print(k, v)

with open(in_file) as handler:
    genes_attr = None

    if file_type == "saf":
        genes_attr = get_saf(handler)
    if file_type == "gtf":
        genes_attr = get_gtf(handler, feat_types)
    if file_type == "gff":
        genes_attr = get_gff(handler, feat_types, biotype_id_mapping=biotype_id_mapping)

    if genes_attr is not None:

        header = True
        for chrom, attr in sorted(genes_attr.items()):
            if header:
                print('\t'.join(("Gene.ID", "Chrom", "Gene.Name", "Biotype")))
                header = False

            for gene_id, values in sorted(attr.items()):
                print('\t'.join((gene_id, chrom, values["gene_name"], values["biotype"])))

    else:
        sys.exit(
            "ERROR: This should have happened. check your gene annotation "
            "and gene_idx.txt files")
