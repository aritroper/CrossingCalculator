from enum import Enum, auto
from itertools import product
import copy

class Gender(Enum):
    MALE = auto()
    FEMALE = auto()

    def __str__(self):
        if self == Gender.MALE:
            return "(M)"
        elif self == Gender.FEMALE:
            return "(F)"
        else:
            return "(M/F)"

class Phenotype(Enum):
    WHITE_EYES = auto()
    ORANGE_EYES = auto()
    RED_EYES = auto()
    STAR_EYES = auto()
    CURLY_WINGS = auto()
    SHORT_BRISTLES = auto()
    BIG_HALTERES = auto()
    SHOULDER_HAIR = auto()
    MALE = auto()
    NONE = auto()

    def __str__(self):
        return ' '.join(word.capitalize() for word in self.name.split('_'))

    # Implement the less than method for sorting
    def __lt__(self, other):
        if isinstance(other, Phenotype):
            return self.value < other.value
        return NotImplemented

class Allele:
    def __init__(self, name, phenotype, recombination, homozygous_lethal, visible):
        self.name = name
        self.phenotype = phenotype
        self.recombination = recombination
        self.homozygous_lethal = homozygous_lethal
        self.visible = visible
        self.recombinated_with = None

    def __str__(self):
        return self.name

    def get_phenotype(self, has_white_eyes):
        return self.phenotype

    def __eq__(self, other):
        if isinstance(other, Allele):
            if self.recombinated_with and other.recombinated_with:
                return sorted([self.name, self.recombinated_with.name]) == sorted([other.name, other.recombinated_with.name])
            elif not self.recombinated_with and not other.recombinated_with:
                return self.name == other.name
            else:
                return False
        return NotImplemented

    def __hash__(self):
        return hash(self.name)

    def recombinate_with(self, a):
        if isinstance(a, Allele):
            self.recombinated_with = a
            return self
        else:
            raise ValueError("Allele recombinate_with: Must recombinate with allele")

class Marker(Allele):
    def __init__(self, name, phenotype):
        super().__init__(name, phenotype, True, True, True)

class Balancer(Allele):
    def __init__(self, name, phenotype):
        super().__init__(name, phenotype, False, True, True)

class Transgene(Allele):
    def __init__(self, name):
        super().__init__(name, Phenotype.NONE, True, False, False)

    def get_phenotype(self, has_white_eyes):
        if has_white_eyes:
            self.visible = True
            return Phenotype.ORANGE_EYES
        else:
            self.visible = False
            return Phenotype.NONE

class Mutation(Allele):
    def __init__(self, name, phenotype):
        super().__init__(name, phenotype, True, False, True)

class AlleleType():

    Y = Allele("Y", Phenotype.MALE, False, True, True)
    PLUS = Allele("PLUS", Phenotype.NONE, True, False, False)

    # Markers
    SM6 = Marker("SM6", Phenotype.CURLY_WINGS)
    S = Marker("S", Phenotype.STAR_EYES)
    UBX = Marker("UBX", Phenotype.BIG_HALTERES)
    TB = Marker("TB", Phenotype.SHOULDER_HAIR)

    # Balancers
    CyO = Balancer("CyO", Phenotype.CURLY_WINGS)
    TM2 = Balancer("TM2", Phenotype.BIG_HALTERES)
    TM3 = Balancer("TM3", Phenotype.SHORT_BRISTLES)
    TM6b = Balancer("Tm6b", Phenotype.SHOULDER_HAIR)

    # Transgenes
    XGAL4 = Transgene("XGAL4")
    DNAD = Transgene("DN-AD")
    SOS = Transgene("SOS")
    DNDB = Transgene("DN-DB")
    ORB = Transgene("ORB")
    UAS = Transgene("UAS-SPARC")

    # Mutations
    WMINUS = Mutation("W-", Phenotype.WHITE_EYES)
    WPLUS = Mutation("W+", Phenotype.RED_EYES)

class Gene():
    def __init__(self, allele1, allele2=None):
        if allele2 is None:
            # Only one gene provided, set gene2 to a default value
            if allele1.homozygous_lethal:
                allele2 = AlleleType.PLUS
            else:
                allele2 = allele1
        if (allele1 == allele2) and (allele1.homozygous_lethal or allele2.homozygous_lethal):
            raise ValueError("Gene Error: allele1 and allele2 must be different")
        if not isinstance(allele1, Allele) or not isinstance(allele2, Allele):
            raise ValueError("Gene Error: allele1 and allele2 must be instances of allele")

        self.allele1 = allele1
        self.allele2 = allele2

    def get_allele_not_transgene(self):
        if isinstance(self.allele1, Transgene) and isinstance(self.allele2, Transgene):
            return None
        elif isinstance(self.allele1, Transgene):
            return self.allele2
        else:
            return self.allele1

    def get_phenotype(self, has_white_eyes):
        if (self.allele1.recombinated_with):
            return sorted([self.allele1.get_phenotype(has_white_eyes), self.allele1.recombinated_with.get_phenotype(has_white_eyes), self.allele2.get_phenotype(has_white_eyes)]) 
        else:
            return sorted([self.allele1.get_phenotype(has_white_eyes), self.allele2.get_phenotype(has_white_eyes)])

    def can_undergo_recombination(self):
        return self.allele1.recombination and self.allele2.recombination and (self.allele1 != self.allele2)

    # Gene undergoes recombination
    def undergo_recombination(self):
        if self.can_undergo_recombination():
            recombinated_allele = copy.deepcopy(self.allele1)
            recombinated_allele.recombinate_with(copy.deepcopy(self.allele2))
            recombinated_allele.recombination = False
            return Gene(recombinated_allele, AlleleType.PLUS)
        else:
            raise ValueError("Gene Error: this gene cannot undergo recombination-- do not force it to!")

    def __str__(self):
        if (self.allele1 == self.allele2):
            return f"{self.allele1.name}"
        elif (self.allele1 == AlleleType.PLUS and self.allele2 != AlleleType.PLUS):
            return f"{self.allele2.name}/{self.allele1.name}"
        elif (self.allele1.recombinated_with != None):
             return f"{self.allele1.name} {self.allele1.recombinated_with.name}/{self.allele2.name}"
        elif (self.allele2 == AlleleType.Y):
            return f"{self.allele1.name}"
        elif (self.allele1 == AlleleType.Y):
            return f"{self.allele2.name}"
        else:
            return f"{self.allele1.name}/{self.allele2.name}"

    def __eq__(self, other):
        return ((self.allele1, self.allele2) == (other.allele1, other.allele2)) or ((self.allele2, self.allele1) == (other.allele1, other.allele2))

    def __hash__(self):
        return hash(self.allele1) + hash(self.allele2)

class Line:
    def __init__(self, gene1=AlleleType.PLUS, gene2=AlleleType.PLUS, gene3=AlleleType.PLUS, gene4=AlleleType.PLUS, id=-1):
        genotype = [gene1, gene2, gene3, gene4]

        for i in range(len(genotype)):
            # Check if the gene is an instance of Allele and has homozygous_lethal=True
            if isinstance(genotype[i], Allele):
                genotype[i] = Gene(genotype[i])

        self.genotype = genotype
        self.has_white_eyes = (self[0] == Gene(AlleleType.WMINUS, AlleleType.WMINUS)) or (self[0] == Gene(AlleleType.WMINUS, AlleleType.Y))
        self.id = id
        
        if (self[0].allele1 == AlleleType.Y or self[0].allele2 == AlleleType.Y):
            self.gender = Gender.MALE
        else:
            self.gender = Gender.FEMALE

    def can_undergo_recombination(self):
        can_undergo_recombination = False
        for i in range(4):
            can_undergo_recombination = can_undergo_recombination or self[i].can_undergo_recombination()
        return can_undergo_recombination and (self.gender == Gender.FEMALE)

    def undergo_recombination(self):
        if self.can_undergo_recombination():
            for i in range(4):
                if (self[i].can_undergo_recombination()):
                    self.genotype[i] = self[i].undergo_recombination()
        else:
            raise ValueError("Line Error: this line cannot undergo recombination-- do not force it to!")

    def get_genotype(self):
        """Return a list of genes in the starter line."""
        return self.genotype

    def get_visual_alleles(self):
        alleles = set()
        for i in range(len(self.genotype)):
            if not isinstance(self[i], Gene):
                raise ValueError("All genes in genotype must be of type gene")
            else:
                if (self[i].allele1.visible):
                    alleles.add(self[i].allele1)
                if (self[i].allele2.visible):
                    alleles.add(self[i].allele2)
        return alleles

    def get_phenotype(self):
        phenotype = []
        for gene in self.genotype:
            phenotype.append(gene.get_phenotype(self.has_white_eyes))
        return phenotype

    def __getitem__(self, index):
        """Get the gene at the specified index."""
        if 0 <= index < len(self.genotype):
            return self.genotype[index]
        else:
            raise IndexError("Index out of range")

    def __eq__(self, other):
        return (self[0] == other[0] and self[1] == other[1] and self[2] == other[2] and self[3] == other[3])

    def __hash__(self):
        return (hash(self[0])) + (hash(self[1])) + (hash(self[2])) + (hash(self[3]))

    def __str__(self):
        """String representation of the starter line."""
        gene_strs = []
        for gene in self.genotype:
            if isinstance(gene, Gene):
                gene_strs.append(str(gene))
            else:
                gene_strs.append("Unknown")  # Fallback for unexpected types

        return ", ".join(gene_strs)


def cross_chromosome(gene1, gene2):
        if not isinstance(gene1, Gene) or not isinstance(gene2, Gene):
            raise ValueError("cross_chromosome Error: all genes should be represented as two alleles")

        crossed = set()  # Initialize an empty set
        alleles1 = [gene1.allele1, gene1.allele2]
        alleles2 = [gene2.allele1, gene2.allele2]

        # Perform initial crossing
        for a1 in alleles1:
            for a2 in alleles2:
                try:
                    g = Gene(a1, a2)
                    crossed.add(g)

                except ValueError:
                    continue

        return crossed

def cross_lines(line1, line2):
     # Return an empty list if the gender condition is not met
    if not ((line1.gender == Gender.MALE and line2.gender == Gender.FEMALE) or (line1.gender == Gender.FEMALE and line2.gender == Gender.MALE)):
        return []

    all_chromosome_gene_combinations = []

    for i in range(4):
        crossed_genes = cross_chromosome(line1[i], line2[i])
        all_chromosome_gene_combinations.append(crossed_genes)

    lines = [Line(*combination) for combination in product(*all_chromosome_gene_combinations)]

    # Handling non-hashable phenotypes
    phenotype_dict = {}
    for line in lines:
        if (line.can_undergo_recombination()):
            line.undergo_recombination()

        phenotype = line.get_phenotype()
        phenotype_key = str(phenotype)  # Convert phenotype to a string or other hashable form

        if phenotype_key in phenotype_dict:
            phenotype_dict[phenotype_key] = None  # Mark phenotype as duplicated
        elif phenotype_key not in phenotype_dict:
            phenotype_dict[phenotype_key] = line

    # Filter out lines with duplicated phenotypes
    unique_lines = [line for line in phenotype_dict.values() if line is not None]

    return unique_lines

'''
Returns False if the target contains alleles not present in its starter lines.
'''
def dummy_check(starter_lines, target):
     # Collect all alleles from the starter lines
    starter_alleles = set()
    for line in starter_lines:
        for gene in line.genotype:
            starter_alleles.add(gene.allele1)
            starter_alleles.add(gene.allele2)

    # Check if all alleles in the target line are present in the starter lines
    for gene in target.genotype:
        if gene.allele1 not in starter_alleles or gene.allele2 not in starter_alleles:
            return False  # Target contains an allele not found in starter lines

    return True

def compute_crosses(starter_lines, target):
    seen_lines = set(starter_lines)
    queue = [(line, []) for line in starter_lines]

    if not (dummy_check(starter_lines, target)):
        return None

    while queue:
        new_queue = []
        added_any_new_line = False  # Flag to check if new lines are added

        for i, (line1, path1) in enumerate(queue):
            for j, (line2, path2) in enumerate(queue):
                if i == j:  # Skip crossing a line with itself
                    continue

                crossed_lines = cross_lines(line1, line2)

                for crossed_line in crossed_lines:
                    if crossed_line in seen_lines:
                        continue  # Skip if line has already been seen

                    if crossed_line == target:
                        return path1 + path2 + [(line1, line2, crossed_line)]

                    new_queue.append((crossed_line, path1 + path2 + [(line1, line2, crossed_line)]))
                    seen_lines.add(crossed_line)
                    added_any_new_line = True

        if not added_any_new_line:
            return None  # No new lines added, target line is unattainable

        queue = new_queue

def compute_pick_for(line1, line2, end_line):
    alleles1 = line1.get_visual_alleles()
    alleles2 = line2.get_visual_alleles()

    all_alleles = alleles1.union(alleles2)

    pick_for = set()
    pick_against = set()

    alleles_end = end_line.get_visual_alleles()

    # Compute pick-for, pick-against
    for allele in all_alleles:
        if allele in alleles_end:
            pick_for.add(allele)
        else:
            pick_against.add(allele)

    # Print pick-for, pick-against
    for a in pick_for:
        print("+" + str(a))

    for a in pick_against:
        print("-" + str(a))

def init_starter_lines(starter_lines):
    LINE_ID = 1
    modified_starter_lines = []

    for line in starter_lines:
        # Copy the original line
        original_line = copy.deepcopy(line)
        original_line.id = LINE_ID

        modified_starter_lines.append(original_line)

        # Create a copy with the opposite gender
        opposite_gender_line = copy.deepcopy(line)

        if line.gender == Gender.MALE:
            # If original is male, create a female copy
            opposite_gender_line.gender = Gender.FEMALE
            # For females, second allele of the first gene should not be Y
            if isinstance(opposite_gender_line[0], Gene) and opposite_gender_line[0].allele2 == AlleleType.Y:
                opposite_gender_line.genotype[0] = Gene(opposite_gender_line[0].allele1, opposite_gender_line[0].allele1)
                if opposite_gender_line.can_undergo_recombination():
                    raise ValueError("Init starter lines-- starting lines should not be able to undergo recombination")
            else:
                raise ValueError("Init starter lines-- creating opposite gender error")
        else:
            # If original is female, create a male copy
            opposite_gender_line.gender = Gender.MALE
            # For males, second allele of the first gene should be Y
            if isinstance(opposite_gender_line[0], Gene):
                opposite_gender_line.genotype[0] = Gene(opposite_gender_line[0].allele1, AlleleType.Y)
            else:
                raise ValueError("Init starter lines-- creating opposite gender error")
        
        opposite_gender_line.id = LINE_ID
        LINE_ID += 1
        modified_starter_lines.append(opposite_gender_line)

    return modified_starter_lines, LINE_ID

def print_crosses(crosses, num_starter_lines):
    if crosses == None:
        print("NO SOLUTION FOUND")
        return

    CROSS = 1
    LINE_ID = num_starter_lines

    print(" ")

    for lines in crosses:

        # Assign these lines ids 
        for i in range(3):
            if (lines[i].id == -1):
                lines[i].id = LINE_ID
                LINE_ID += 1

        print("CROSS " + str(CROSS))
        print("(" + str(lines[0].id) + ") " + str(lines[0]) + " | " + str(lines[0].gender))
        print(" X")
        print("(" + str(lines[1].id) + ") " + str(lines[1]) + " | " + str(lines[1].gender))
        print(" ---------------------")
        print("(" + str(lines[2].id) + ") " + str(lines[2]))
        print(" ")

        compute_pick_for(lines[0], lines[1], lines[2])

        print(" ")

        CROSS += 1

def main(starter_lines, target):
    s_lines = init_starter_lines(starter_lines)
    path = compute_crosses(s_lines[0], target)
    print_crosses(path, s_lines[1])

# Starter lines
lineA = Line(AlleleType.WPLUS, AlleleType.PLUS, Gene(AlleleType.TM2, AlleleType.TM6b), AlleleType.PLUS)
lineB = Line(Gene(AlleleType.WPLUS, AlleleType.Y), Gene(AlleleType.S, AlleleType.CyO), Gene(AlleleType.XGAL4, AlleleType.TM6b), AlleleType.PLUS)
lineC = Line(AlleleType.WPLUS, AlleleType.PLUS, Gene(AlleleType.DNAD, AlleleType.TM6b), AlleleType.PLUS)

target = Line(AlleleType.WPLUS, AlleleType.PLUS, Gene((AlleleType.DNAD).recombinate_with(AlleleType.XGAL4), AlleleType.TM6b), AlleleType.PLUS)

main([lineA, lineB, lineC], target)
