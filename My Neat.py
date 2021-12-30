import random
import math
import copy

class Node:
    def __init__(self, number, position):
        self.number = number
        self.position = position
        self.value = 0.0

class Gene:
    def __init__(self, innode, outnode):
        self.innode = innode
        self.outnode = outnode
        self.active = True
        self.weight = random.uniform(-0.5, 0.5)
        innum = innode.number
        outnum = outnode.number
        self.innovation = (0.5 * ((innum + outnum)*(innum + outnum+1))) + outnum  #Cantor pairing function
    def to_string(self):
        print("Inno: ", self.innovation, " Innode: ", self.innode.number, " Outnode: ", self.outnode.number, " Weight: ", self.weight, " Active: ", self.active)

class Genome:
    def __init__(self, num_input, num_output, create_genes=True):
        self.num_inputs = num_input
        self.num_outputs = num_output
        self.fitness = 0
        self.alive = True
        self.node_count = int(1)
        self.gene_list = []

        if create_genes:
            in_list = []
            for _ in range(num_input + 1):
                in_list.append(Node(self.node_count, 0 ))
                self.node_count += 1

            out_list = []
            for _ in range(num_output):
                out_list.append(Node(self.node_count, 4096))
                self.node_count += 1

            for i in in_list:
                for o in out_list:
                    self.gene_list.append(Gene(i, o))

            del in_list
            del out_list
    
    def crossover(self, genome):
        baby = Genome(self.num_inputs, self.num_outputs, create_genes=False)
        gene_copy = copy.deepcopy(self.gene_list)
        baby.node_count = self.node_count
        
        this_inno = [gene.innovation for gene in self.gene_list]
        that_inno = [gene.innovation for gene in genome.gene_list]

        for i, num in enumerate(this_inno):
            gene = gene_copy[i]
            if (num in that_inno) and (random.random() > 0.5):
                index = that_inno.index(num)
                gene.weight = genome.gene_list[index].weight
                gene.active = genome.gene_list[index].active
                if random.random() > 0.75:
                    gene.active = True
            baby.gene_list.append(gene)
                
        assert len(baby.gene_list) == len(self.gene_list), "Crossover is not working properly"
        return baby

    def mutate_node(self):
        old_gene =  random.choice(self.gene_list)
        old_gene.active = False
        position = (old_gene.innode.position + old_gene.outnode.position)/2
        new_node = Node(self.node_count, position)
        self.node_count += 1

        gene_1 = Gene(old_gene.innode, new_node)
        gene_1.weight = old_gene.weight

        gene_2 = Gene(new_node, old_gene.outnode)
        gene_2.weight = 1

        self.gene_list.append(gene_1)
        self.gene_list.append(gene_2)

        del gene_1
        del gene_2
        del old_gene

    def mutate_connection(self):
        zipped = [(gene.innode, gene.outnode) for gene in self.gene_list]

        valid = False
        innode, outnode = None, None
        for _ in range(20):
            innode = zipped[random.randint(0, len(zipped)-1)][0]
            outnode = zipped[random.randint(0, len(zipped)-1)][1]
            if outnode.position < innode.position:
                innode, outnode = outnode, innode
            if ((innode, outnode) not in zipped) and (outnode.position != innode.position):
                valid = True
                break
        if valid:
            new_gene = Gene(innode, outnode)
            self.gene_list.append(new_gene)
        
        del zipped
        del innode
        del outnode

    def mutate_weights(self):
        for g in self.gene_list:
            if random.random() < 0.9:
                g.weight += random.uniform(-1, 1) *0.06
            else:
                g.weight = random.random() * 2 -1
            if g.weight > 1:
                g.weight = 1
            elif g.weight < -1:
                g.weight = -1      

    def mutate_active(self):
        for g in self.gene_list:
            if random.random() < 0.25:
                g.active = not g.active
    
    def feed_forward(self, inputs):
        assert len(inputs) == self.num_inputs, "Not enough arguements for the genome!"
        output_nodes = []
        input_nodes = []
        outputs = []
        x = [1]
        x.extend(list(inputs))

        self.gene_list = sorted(self.gene_list, key=lambda Gene: Gene.innode.number)
        for gene in self.gene_list:
            if gene.innode.position == 0 and gene.innode not in input_nodes:
                input_nodes.append(gene.innode)
            if len(input_nodes) == self.num_inputs+1:
                break

        for i, node in enumerate(input_nodes):
            node.value = x[i]

        network = [gene for gene in self.gene_list if gene.active]
        network = sorted(network, key=lambda Gene: Gene.innode.position)
        for gene in network:
            inn = gene.innode
            out = gene.outnode
            if inn.position == 0:
                out.value += inn.value * gene.weight
            else:
                out.value += self.activation(inn.value * gene.weight)

        self.gene_list = sorted(self.gene_list, key=lambda Gene: Gene.outnode.position, reverse=True)
        max_pos = self.gene_list[0].outnode.position
        for gene in self.gene_list:
            if gene.outnode.position < max_pos:
                gene.outnode.value = 0
            elif gene.outnode not in output_nodes:
                output_nodes.append(gene.outnode)
                outputs.append(self.activation(gene.outnode.value))
                gene.outnode.value = 0
                
        assert len(outputs) == self.num_outputs, "Genome output is incorrect"
        return outputs
 
    def activation(self, x):
        denom = (1+math.exp(-4.9*x))
        return 1/denom    
        
class Species:
    def __init__(self, genome, disjoint_coef, excess_coef, weight_coef, speciation_threshold, stagnant_buffer):
        self.genomes = [genome]
        self.best = genome
        self.rep = genome
        self.avg_fitness = 0
        self.prev_fitness = 0
        self.stagnation = int(0)
        self.stagnant_buffer = stagnant_buffer
        self.spawn_num = int(0)

        self.excess_coef = excess_coef
        self.disjoint_coef = disjoint_coef
        self.weight_coef = weight_coef
        self.speciation_threshold = speciation_threshold

    def speciate(self, genome_1):
        excess = int(0)
        disjoint = int(0)
        matching = int(0)
        weight_diff = 0

        g1_inno = [gene.innovation for gene in genome_1.gene_list]
        g2_inno = [gene.innovation for gene in self.rep.gene_list]

        for i, inno_1 in enumerate(g1_inno):
            for j, inno_2 in enumerate(g2_inno):
                if inno_1 == inno_2:
                    matching += 1
                    weight_diff += abs(genome_1.gene_list[i].weight - self.rep.gene_list[j].weight)
                    break
        
        highest = max(g2_inno)
        for inno in g1_inno:
            if (inno not in g2_inno) and (inno < highest):
                disjoint +=1
            elif (inno not in g2_inno) and (inno > highest):
                excess += 1

        N = max(len(g1_inno), len(g2_inno))
        if N < 20: N = 1
        threshold = (self.excess_coef*excess/N) + (self.disjoint_coef*disjoint/N) + ((weight_diff/matching)*self.weight_coef)
        #del g1_inno, g2_inno, excess, matching, disjoint, weight_diff
        if threshold < self.speciation_threshold: 
            self.genomes.append(genome_1)
            return True
        else: return False

    def calculate_spawn(self, avg_fit, pop_count):
        self.spawn_num = int(self.avg_fitness / avg_fit * pop_count)-1
    
    def make_baby(self):
        baby = None
        if random.random() > 0.8:
            baby = copy.deepcopy(self.select_genome())
        else:
            parent_1 = self.genomes[0]
            parent_2 = self.select_genome()
            if parent_2.fitness > parent_1.fitness:
                baby = parent_2.crossover(parent_1)
            else: baby = parent_1.crossover(parent_2)
        return baby
    
    def select_genome(self):
        fitness_sum = sum([g.fitness for g in self.genomes])
        rando = random.randint(0, int(fitness_sum))
        running_sum = 0

        for g in self.genomes:
            running_sum += g.fitness
            if running_sum > rando:
                return g
        else:
            return self.genomes[0]

    def set_avg_fitness(self):
        self.avg_fitness = sum([g.fitness for g in self.genomes]) / len(self.genomes)

    def cull_genomes(self, cull_rate):
        self.genomes = sorted(self.genomes, key=lambda Genome: Genome.fitness, reverse=True)
        length = len(self.genomes)
        num_dead = int(len(self.genomes) * cull_rate)
        self.genomes = self.genomes[:(length - num_dead)]
        assert len(self.genomes) > 0, "Killing all genomes with culling: " + str(num_dead) + " "+ str(length) 

    def fitness_sharing(self):
        for g in self.genomes:
            g.fitness /= len(self.genomes)

    def check_stagnation(self):
        if self.best.fitness > self.prev_fitness + self.stagnant_buffer:
            self.stagnation = 0
            self.prev_fitness = self.best.fitness
        else:
            self.stagnation += 1

    def go_extinct(self):
        self.stagnantion = -99
        self.genomes.clear()
        self.best = None
        self.rep = None
        self.fitness = 0

class Population:
    def __init__(self, count, num_inputs, num_outputs):
        self.generation = int(1)
        self.players = []
        self.count = count
        self.num_inputs = num_inputs
        self.num_outputs = num_outputs
        self.mass_extinction = False
        
        self.species = []
        self.champion = None
        self.stagnant_max = int(20)
        self.stagnant_buffer = 0.08
        self.cull_rate = 0.6
        self.excess = 1
        self.disjoint = 1
        self.weight = 0.5
        self.thresh = 3

        self.weight_rate = 0.8
        self.node_rate = 0.03
        self.gene_rate = 0.5
        self.enable_rate = 0.25

        for _ in range(count):
            genome = Genome(self.num_inputs, self.num_outputs)
            self.mutate(genome)
            self.players.append(genome)
        self.champion = self.players[0]
    
    def make_list(self):
        assert len(self.players) == self.count, "Actual: "+ str(len(self.players)) + " Wanted: " + str(self.count)
        return self.players

    def mutate(self, genome):
        if random.random() < self.weight_rate: genome.mutate_weights()
        if random.random() < self.node_rate: genome.mutate_node()
        if random.random() < self.gene_rate: genome.mutate_connection()
        if random.random() < self.enable_rate: genome.mutate_active()

    def evolve(self):
        self.place_genomes()
        self.sort_and_cull_species()
        self.set_champion()
        self.kill_species()
        if self.mass_extinction:
            self.repopulate()
            return False

        children = []
        avg_fit_sum = sum([s.avg_fitness for s in self.species])
        for s in self.species:
            children.append(copy.deepcopy(s.best))
            no_of_children = int(s.avg_fitness / avg_fit_sum * len(self.players)) - 1
            for _ in range(no_of_children):
                child = s.make_baby()
                self.mutate(child)
                children.append(child)
        for _ in range(self.count - len(children)):
            child = self.species[0].make_baby()
            self.mutate(child)
            children.append(child)

        self.players.clear()
        self.generation += 1
        self.players = [child for child in children]
        del children
        return True
            
    def place_genomes(self):
        for s in self.species:
            s.rep = s.select_genome()
            s.genomes.clear()
        for genome in self.players:
            for s in self.species:
                placed = s.speciate(genome)
                if placed: break
            else:
                specie = Species(genome, self.disjoint, self.excess, self.weight, self.thresh, self.stagnant_buffer)
                self.species.append(specie)
        i = 0
        for _ in range(len(self.species)):
            s = self.species[i]
            if not s.genomes: self.species.pop(i)
            else: i += 1

        for s in self.species:
            assert s.genomes, "A species with no genomes made it through"

    def sort_and_cull_species(self):
        for s in self.species:
            s.genomes = sorted(s.genomes, key=lambda Genome: Genome.fitness, reverse=True)
            s.cull_genomes(self.cull_rate)
            s.best = s.genomes[0]
        self.species = sorted(self.species, key=lambda Species: Species.best.fitness, reverse=True)
        
    def set_champion(self):
        for s in self.species:
            if s.best.fitness > self.champion.fitness:
                self.champion = s.best              

    def kill_species(self):
        avg_fit_sum = 0
        for s in self.species:
            s.set_avg_fitness()
            avg_fit_sum += s.avg_fitness

        index = 0
        for _ in range(len(self.species)):
            s = self.species[index]
            s.check_stagnation()
            if s.stagnation > self.stagnant_max:
                s.go_extinct()
                self.species.pop(index)
            elif (s.avg_fitness / avg_fit_sum * len(self.players)) - 1 < 1:
                s.go_extinct()
                self.species.pop(index)
            else: index += 1

        if not self.species:
            self.mass_extinction = True   

    def repopulate(self):
        #print("Mass extinction at generation: ", self.generation)
        #print("Repopulating with fresh genomes")
        self.players.clear()
        for _ in range(self.count):
            genome = Genome(self.num_inputs, self.num_outputs)
            self.mutate(genome)
            self.players.append(genome)
        self.champion = self.players[0]
        self.mass_extinction = False




#Kenneth Stanley said this on xor: "I generally said < 0.5 is 0 and >= 0.5 is 1 for the purposes of a solution"
#I will be using this to check if a solution has been found

#TODO The avg fitness of genomes is not improving, need to find possible causes
#IT'S DONE! It occured twice in a row! The algorithm still isnt perfect and requires multiple attempts. Genomes also have too many genes (node_rate 0.03, gene_rate 0.05)
#Changed crossover so that parent one is always the fittest genome and select player is back to normal
#Allowed for the most efficient solution possible, 7 genes in 27 generations!

#Genomes are still mutating too rapidly, need to change species offspring to some form of distribution relative to their number and fitness
#Changed stagnantion to go by avg_fitness
#AVG's improved as the cull rate was set to 0.2

#Out of 3 attempts of 50, completed was 16, 18, 18
#changed cull rate to 0.8: completed was 13, 17, 17
#Mass extinctions still happening lots of times. A full run usually either ends in extinction or finding.
#Will come up with an algorithm for more evnly splitting offspring between species


#pop = Population(200,2,1)

completed = 0
uncompleted = 0
avg_found_len = 0
avg_generations_found = 0
mass_extinctions = 0

for _ in range(30):
    pop = Population(200,2,1)
    found = False
    i = 0
    while i < 100 and not found:
        i += 1
        #print("-------Generation: ",i,"---------")
        gl = pop.make_list()
        for g in pop.players:
            ones = []
            zeros = []

            out = g.feed_forward([1,1])
            ii =  abs(0 - out[0])
            zeros.append(out[0])

            out = g.feed_forward([1,0])
            io = abs(1 - out[0])
            ones.append(out[0])
            
            out = g.feed_forward([0,0])
            oo = abs(0 - out[0])
            zeros.append(out[0])

            out = g.feed_forward([0,1])
            oi = abs(1 - out[0])
            ones.append(out[0])

            if zeros[0] < 0.5 and zeros[1] < 0.5 and ones[0] >= 0.5 and ones[1] >= 0.5: 
                found = True
                avg_found_len += len(g.gene_list)
                #print("\nA solution has been found! ",zeros, ones)
                #for gene in g.gene_list:
                #   gene.to_string()
                break
            score = 4 - (ii + oo + io + oi)
            score *= score
            g.fitness = score
        if found:
            #print("A genome has successfully passed the test")
            completed += 1
            avg_generations_found += i
        else:
            not_extinct = pop.evolve()
            if not not_extinct:
                mass_extinctions += 1
                break
    if not found:
        uncompleted += 1
    del pop.players
    del pop.species
    del pop

print("Uncompleted: ", uncompleted)
print("Mass extinctions: ", mass_extinctions)
print("Completed: ", completed)
print("Avg completion length: ", avg_found_len/completed)
print("Avg gens to completion: ", avg_generations_found/completed)


#best = max([g.fitness for g in gl])
#print("Highest scoring genome: ", best)
#print("Species avg fitnesses: ", [s.avg_fitness for s in pop.species])    
#print()





