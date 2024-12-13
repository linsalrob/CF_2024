
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

    def get_taxid(self, pathogen):
        """
        Get the taxonomic ID of a bacterial pathogen.
        :param pathogen: str
        :return: int
        """
        return self.pathogens.get(pathogen, None)