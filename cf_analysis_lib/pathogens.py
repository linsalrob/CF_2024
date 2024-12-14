
"""

"""


class BacterialPathogens:
    """
    Class to store bacterial pathogens.
    """

    def __init__(self):
        self.pathogens = {
            "Acinetobacter baumannii": 470,
            "Actinomyces israelii": 1659,
            "Bacillus anthracis": 1392,
            "Bordetella pertussis": 520,
            "Borrelia burgdorferi": 138,
            "Brucella melitensis": 29459,
            "Burkholderia cenocepacia": 95486,
            "Campylobacter jejuni": 197,
            "Chlamydia trachomatis": 813,
            "Clostridium botulinum": 1491,
            "Clostridium difficile": 1497,
            "Clostridium perfringens": 1502,
            "Corynebacterium diphtheriae": 1717,
            "Enterococcus faecalis": 1350,
            "Escherichia coli": 562,
            "Francisella tularensis": 263,
            "Haemophilus influenzae": 727,
            "Helicobacter pylori": 210,
            "Klebsiella pneumoniae": 573,
            "Legionella pneumophila": 446,
            "Listeria monocytogenes": 1639,
            "Moraxella": 475,
            "Mycobacterium tuberculosis": 1773,
            "Mycoplasma": 2093,
            "Neisseria gonorrhoeae": 485,
            "Neisseria meningitidis": 487,
            "Nocardia abscessus":  120957,
            "Pasteurella multocida" : 747,
            "Pseudomonas aeruginosa": 287,
            "Rickettsia prowazekii": 782,
            "Salmonella enterica": 289,
            "Shigella dysenteriae": 622,
            "Staphylococcus aureus": 1280,
            "Streptococcus pneumoniae": 1313,
            "Streptococcus pyogenes": 1314,
            "Vibrio cholerae": 666,
            "Yersinia pestis": 632,
        }



        self.taxonomy_strings = {
            '487' : 'k__Bacteria;p__Pseudomonadota;c__Betaproteobacteria;o__Neisseriales;f__Neisseriaceae;g__Neisseria;s__Neisseria meningitidis',
            '1392' : 'k__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Bacillaceae;g__Bacillus;s__Bacillus anthracis',
            '446' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Legionellales;f__Legionellaceae;g__Legionella;s__Legionella pneumophila',
            '1314' : 'k__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus pyogenes',
            '727' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Haemophilus;s__Haemophilus influenzae',
            '480' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Moraxellales;f__Moraxellaceae;g__Moraxella;s__Moraxella catarrhalis',
            '1351' : 'k__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Enterococcaceae;g__Enterococcus;s__Enterococcus faecalis',
            '1313' : 'k__Bacteria;p__Bacillota;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus pneumoniae',
            '29459' : 'k__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Hyphomicrobiales;f__Brucellaceae;g__Brucella;s__Brucella melitensis',
            '782' : 'k__Bacteria;p__Pseudomonadota;c__Alphaproteobacteria;o__Rickettsiales;f__Rickettsiaceae;g__Rickettsia;s__Rickettsia prowazekii',
            '1491' : 'k__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Clostridium;s__Clostridium botulinum',
            '813' : 'k__Bacteria;p__Chlamydiota;c__Chlamydiia;o__Chlamydiales;f__Chlamydiaceae;g__Chlamydia;s__Chlamydia trachomatis',
            '622' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Shigella;s__Shigella dysenteriae',
            '139' : 'k__Bacteria;p__Spirochaetota;c__Spirochaetia;o__Spirochaetales;f__Borreliaceae;g__Borreliella;s__Borreliella burgdorferi',
            '1659' : 'k__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Actinomycetales;f__Actinomycetaceae;g__Actinomyces;s__Actinomyces israelii',
            '1717' : 'k__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales;f__Corynebacteriaceae;g__Corynebacterium;s__Corynebacterium diphtheriae',
            '1280' : 'k__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus',
            '210' : 'k__Bacteria;p__Campylobacterota;c__Epsilonproteobacteria;o__Campylobacterales;f__Helicobacteraceae;g__Helicobacter;s__Helicobacter pylori',
            '1497' : 'k__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Clostridium;s__Clostridium formicaceticum',
            '289' : 'k__Bacteria;p__Pseudomonadota;c__Betaproteobacteria;o__Burkholderiales;f__Comamonadaceae;g__Hydrogenophaga;s__Hydrogenophaga pseudoflava',
            '287' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas aeruginosa',
            '95486' : 'k__Bacteria;p__Pseudomonadota;c__Betaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__Burkholderia cenocepacia',
            '520' : 'k__Bacteria;p__Pseudomonadota;c__Betaproteobacteria;o__Burkholderiales;f__Alcaligenaceae;g__Bordetella;s__Bordetella pertussis',
            '263' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Thiotrichales;f__Francisellaceae;g__Francisella;s__Francisella tularensis',
            '120957' : 'k__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales;f__Nocardiaceae;g__Nocardia;s__Nocardia abscessus',
            '1502' : 'k__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Clostridiaceae;g__Clostridium;s__Clostridium perfringens',
            '1639' : 'k__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;f__Listeriaceae;g__Listeria;s__Listeria monocytogenes',
            '1773' : 'k__Bacteria;p__Actinomycetota;c__Actinomycetes;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__Mycobacterium tuberculosis',
            '562' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli',
            '485' : 'k__Bacteria;p__Pseudomonadota;c__Betaproteobacteria;o__Neisseriales;f__Neisseriaceae;g__Neisseria;s__Neisseria gonorrhoeae',
            '747' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Pasteurellales;f__Pasteurellaceae;g__Pasteurella;s__Pasteurella multocida',
            '2104' : 'k__Bacteria;p__Mycoplasmatota;c__;o__Mycoplasmoidales;f__Mycoplasmoidaceae;g__Mycoplasmoides;s__Mycoplasmoides pneumoniae',
            '470' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Moraxellales;f__Moraxellaceae;g__Acinetobacter;s__Acinetobacter baumannii',
            '573' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Klebsiella;s__Klebsiella pneumoniae',
            '632' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Yersiniaceae;g__Yersinia;s__Yersinia pestis',
            '197' : 'k__Bacteria;p__Campylobacterota;c__Epsilonproteobacteria;o__Campylobacterales;f__Campylobacteraceae;g__Campylobacter;s__Campylobacter jejuni',
            '666' : 'k__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Vibrionales;f__Vibrionaceae;g__Vibrio;s__Vibrio cholerae'
        }

        self.taxonomic_levels = {
            'kingdom' : ['Bacteria'],
            'phylum' : ['Chlamydiota', 'Actinomycetota', 'Pseudomonadota', 'Mycoplasmatota', 'Spirochaetota', 'Bacillota', 'Campylobacterota'],
            'class' : ['Chlamydiia', 'Spirochaetia', 'Epsilonproteobacteria', 'Bacilli', '', 'Betaproteobacteria', 'Gammaproteobacteria', 'Actinomycetes', 'Alphaproteobacteria', 'Clostridia'],
            'order' : ['Vibrionales', 'Moraxellales', 'Bacillales', 'Rickettsiales', 'Campylobacterales', 'Mycoplasmoidales', 'Enterobacterales', 'Actinomycetales', 'Eubacteriales', 'Neisseriales', 'Lactobacillales', 'Spirochaetales', 'Pasteurellales', 'Hyphomicrobiales', 'Mycobacteriales', 'Pseudomonadales', 'Burkholderiales', 'Thiotrichales', 'Chlamydiales', 'Legionellales'],
            'family' : ['Neisseriaceae', 'Corynebacteriaceae', 'Campylobacteraceae', 'Chlamydiaceae', 'Actinomycetaceae', 'Streptococcaceae', 'Rickettsiaceae', 'Nocardiaceae', 'Vibrionaceae', 'Moraxellaceae', 'Francisellaceae', 'Yersiniaceae', 'Comamonadaceae', 'Mycoplasmoidaceae', 'Borreliaceae', 'Legionellaceae', 'Pseudomonadaceae', 'Clostridiaceae', 'Burkholderiaceae', 'Pasteurellaceae', 'Helicobacteraceae', 'Listeriaceae', 'Alcaligenaceae', 'Bacillaceae', 'Staphylococcaceae', 'Brucellaceae', 'Mycobacteriaceae', 'Enterococcaceae', 'Enterobacteriaceae'],
            'genus' : ['Burkholderia', 'Streptococcus', 'Acinetobacter', 'Pasteurella', 'Klebsiella', 'Staphylococcus', 'Listeria', 'Vibrio', 'Francisella', 'Bacillus', 'Moraxella', 'Shigella', 'Helicobacter', 'Pseudomonas', 'Rickettsia', 'Neisseria', 'Mycobacterium', 'Clostridium', 'Borreliella', 'Mycoplasmoides', 'Yersinia', 'Haemophilus', 'Enterococcus', 'Actinomyces', 'Legionella', 'Chlamydia', 'Corynebacterium', 'Escherichia', 'Brucella', 'Hydrogenophaga', 'Campylobacter', 'Nocardia', 'Bordetella'],
            'species' : ['Bacillus anthracis', 'Mycobacterium tuberculosis', 'Nocardia abscessus', 'Burkholderia cenocepacia', 'Brucella melitensis', 'Enterococcus faecalis', 'Streptococcus pyogenes', 'Borreliella burgdorferi', 'Legionella pneumophila', 'Pasteurella multocida', 'Escherichia coli', 'Neisseria gonorrhoeae', 'Listeria monocytogenes', 'Shigella dysenteriae', 'Rickettsia prowazekii', 'Staphylococcus aureus', 'Chlamydia trachomatis', 'Helicobacter pylori', 'Streptococcus pneumoniae', 'Neisseria meningitidis', 'Bordetella pertussis', 'Acinetobacter baumannii', 'Mycoplasmoides pneumoniae', 'Vibrio cholerae', 'Pseudomonas aeruginosa', 'Actinomyces israelii', 'Francisella tularensis', 'Corynebacterium diphtheriae', 'Hydrogenophaga pseudoflava', 'Clostridium formicaceticum', 'Clostridium botulinum', 'Yersinia pestis', 'Campylobacter jejuni', 'Haemophilus influenzae', 'Klebsiella pneumoniae', 'Clostridium perfringens', 'Moraxella catarrhalis']
        }

    
        self.pseudomonas_taxa = {
            'kingdom' : 'Bacteria',
            'phylum' : 'Pseudomonadota',
            'class' : 'Gammaproteobacteria',
            'order' : 'Pseudomonadales',
            'family' : 'Pseudomonadaceae',
            'genus' : 'Pseudomonas',
            'species' : 'Pseudomonas aeruginosa'
        }

    def get_taxid(self, pathogen):
        """
        Get the taxonomic ID of a bacterial pathogen.
        :param pathogen: str
        :return: int
        """
        return self.pathogens.get(pathogen, None)

    def get_taxonomy_string(self, taxid):
        """
        Get the taxonomy string of a bacterial pathogen.
        :param taxid: int
        :return: str
        """
        return self.taxonomy_strings.get(str(taxid), None)

    def get_taxonomic_levels(self, level):
        """
        Get the taxonomic levels of a bacterial pathogen.
        :return: array
        """
        return self.taxonomic_levels.get(level, None)

    def get_pseudomonas_taxa(self, level):
        """
        Get the taxonomic levels of Pseudomonas.
        :return: str
        """
        return self.pseudomonas_taxa.get(level, None)

