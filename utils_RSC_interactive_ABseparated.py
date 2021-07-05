import networkx as nx
from itertools import combinations
import numpy as np
import random
import matplotlib.pyplot as plt


class SimplagionModel():
    def __init__(self, node_neighbors_dict, triangles_list, I_A_percentage, I_B_percentage, I_AB_percentage):
        
        #parameters
        self.neighbors_dict = node_neighbors_dict
        self.triangles_list = triangles_list
        self.nodes = list(node_neighbors_dict.keys())
        self.N = len(node_neighbors_dict.keys())
        self.I_A = I_A_percentage * self.N//100
        self.I_B = I_B_percentage * self.N//100
        self.I_AB = I_AB_percentage * self.N//100
        
        
        #Initial setup
        #I save the infected nodes for each state of the first initialisation in case I want to repeat several runs with
        #the same configuration
        self.infected_this_setup_withA, self.infected_this_setup_withB, self.infected_this_setup_withAB = self.initial_setup()
        
    def initial_setup(self, fixed_nodes_to_infect_withA=None, fixed_nodes_to_infect_withB=None, fixed_nodes_to_infect_withAB=None):
        #going to use this to store the agents in each state
        self.sAgentSet = set()
        self.i_A_AgentSet = set()
        self.i_B_AgentSet = set()
        self.i_AB_AgentSet = set()
        
        #and here we're going to store the counts of how many agents are in each
        #state @ each time step
        self.i_A_List = []
        self.i_B_List = []
        self.i_AB_List = []
        self.S_List = []
        
        self.t = 0
        
        #start with everyone susceptible
        for n in self.nodes:
            self.sAgentSet.add(n)
        
        infected_this_setup_withA=[]
        infected_this_setup_withB=[]
        infected_this_setup_withAB=[]
        #infect nodes
        if fixed_nodes_to_infect_withA==None and fixed_nodes_to_infect_withB==None and fixed_nodes_to_infect_withAB==None : #the first time I create the model (the instance __init__)
            for ite in range(self.I_A): #randomly choosing agents to initially be in state A
                to_infect = random.choice(list(self.sAgentSet))
                self.infectAgent_withA(to_infect)
                infected_this_setup_withA.append(to_infect)
            for ite in range(self.I_B): #randomly choosing agents to initially be in state B
                to_infect = random.choice(list(self.sAgentSet))
                self.infectAgent_withB(to_infect)
                infected_this_setup_withB.append(to_infect)
            for ite in range(self.I_AB): #randomly choosing agents to initially be in state AB
                to_infect = random.choice(list(self.sAgentSet))
                self.infectAgent_withA(to_infect)
                self.infectAgent_withB(to_infect)
                infected_this_setup_withAB.append(to_infect)                
        else: #I already have run the model and this is not the first run, I want to infect the same nodes
            for to_infect in fixed_nodes_to_infect_withA:
                self.infectAgent_withA(to_infect)
                infected_this_setup_withA.append(to_infect)
            for to_infect in fixed_nodes_to_infect_withB:
                self.infectAgent_withB(to_infect)
                infected_this_setup_withB.append(to_infect)
            for to_infect in fixed_nodes_to_infect_withAB:
                self.infectAgent_withA(to_infect)
                self.infectAgent_withB(to_infect)
                infected_this_setup_withAB.append(to_infect)
        self.i_A_List.append(len(self.i_A_AgentSet))
        self.i_B_List.append(len(self.i_B_AgentSet))
        self.i_AB_List.append(len(self.i_AB_AgentSet))
        self.S_List.append(len(self.sAgentSet))
        return infected_this_setup_withA, infected_this_setup_withB, infected_this_setup_withAB 
    
    def infectAgent_withA(self, agent):
        if (agent in self.i_B_AgentSet):
            self.i_B_AgentSet.remove(agent)
            self.i_AB_AgentSet.add(agent)
        if agent in self.sAgentSet:
            self.i_A_AgentSet.add(agent)
            self.sAgentSet.remove(agent)
        return 1
    
    def infectAgent_withB(self, agent):
        if (agent in self.i_A_AgentSet):
            self.i_A_AgentSet.remove(agent)
            self.i_AB_AgentSet.add(agent)
        if agent in self.sAgentSet:
            self.i_B_AgentSet.add(agent)
            self.sAgentSet.remove(agent)
        return 1
    
    def recoverAgent_fromA(self,agent):
        if agent in self.i_AB_AgentSet:
            self.i_AB_AgentSet.remove(agent)
            self.i_B_AgentSet.add(agent)
        if (agent in self.i_A_AgentSet):
            self.i_A_AgentSet.remove(agent)
            self.sAgentSet.add(agent)
        return -1
    
    def recoverAgent_fromB(self,agent):
        if (agent in self.i_AB_AgentSet):
            self.i_AB_AgentSet.remove(agent)
            self.i_A_AgentSet.add(agent)
        if (agent in self.i_B_AgentSet):
            self.i_B_AgentSet.remove(agent)
            self.sAgentSet.add(agent)
        return -1
    
    def run(self, t_max, beta1_A, beta2_A, mu_A, beta1_B, beta2_B, mu_B, epsilon_A, epsilon_B):
        self.t_max = t_max
        while self.t<=self.t_max:
            newI_withA_list = set()
            newI_withB_list = set()
            
            #CONTAGION FOR PATHOGEN A 
            #STANDARD CONTAGION from iAgent in state A
            #we need to loop over the agents who are currently in the state A
            for iAgent in self.i_A_AgentSet:
                #expose their network neighbors
                for agent in self.neighbors_dict[iAgent]:
                    if agent in self.sAgentSet:
                        if (random.random() <= beta1_A):
                            newI_withA_list.add(agent)
                    if (agent in self.i_B_AgentSet):
                        if(random.random() <= beta1_A*epsilon_B):
                            newI_withA_list.add(agent)
            #STANDARD CONTAGION from iAgent in state AB
            #we need to loop over the agents who are currently in the state AB
            for iAgent in self.i_AB_AgentSet:
                #expose their network neighbors
                for agent in self.neighbors_dict[iAgent]:
                    if agent in self.sAgentSet:
                        if (random.random() <= beta1_A):
                            newI_withA_list.add(agent)
                    if (agent in self.i_B_AgentSet):
                        if(random.random() <= beta1_A*epsilon_B):
                            newI_withA_list.add(agent)              
            #TRIANGLE CONTAGION
            for triangle in self.triangles_list:
                n1, n2, n3 = triangle
                if (n1 in self.i_A_AgentSet) or (n1 in self.i_AB_AgentSet):
                    if (n2 in self.i_A_AgentSet) or (n2 in self.i_AB_AgentSet):
                        if n3 in self.sAgentSet: #n1, n2 infected (A or AB) n3 suscectible
                            if (random.random() <= beta2_A):
                                newI_withA_list.add(n3)
                        if (n3 in self.i_B_AgentSet):
                            if(random.random()<= beta2_A*epsilon_B):
                                newI_withA_list.add(agent) #n1, n2 infected n3 in state B
                    else: #n2 in state B or suscectible
                        if ((n3 in self.i_A_AgentSet) or (n3 in self.i_AB_AgentSet)) and (n2 in self.sAgentSet):
                            if (random.random() <= beta2_A): #n1, n3 infected (A or AB) n2 suscectible
                                newI_withA_list.add(n2)
                        if ((n3 in self.i_A_AgentSet) or (n3 in self.i_AB_AgentSet)) and (n2 in self.i_B_AgentSet):
                            if (random.random() <= beta2_A*epsilon_B): #n1, n3 infected (A or AB) n2 in state B
                                newI_withA_list.add(n2)
                else: #n1 in state B or suscectible
                    if ((n2 in self.i_A_AgentSet) or (n2 in self.i_AB_AgentSet)) and ((n3 in self.i_A_AgentSet) or (n3 in self.i_AB_AgentSet)) and (n1 in self.sAgentSet):
                        #n2, n3 infected (A or AB) and n1 suscectible
                        if (random.random() <= beta2_A):
                            newI_withA_list.add(n1)
                    if ((n2 in self.i_A_AgentSet) or (n2 in self.i_AB_AgentSet)) and ((n3 in self.i_A_AgentSet) or (n3 in self.i_AB_AgentSet)) and (n1 in self.i_B_AgentSet):
                        #n2, n3 infected (A or AB) and n1 in state B
                        if (random.random() <= beta2_A*epsilon_B):
                            newI_withA_list.add(n1)
                
            #CONTAGION FOR PATHOGEN B 
            #STANDARD CONTAGION from iAgent in B
            #we need to loop over the agents who are currently in the state B
            for iAgent in self.i_B_AgentSet:
                #expose their network neighbors
                for agent in self.neighbors_dict[iAgent]:
                    if agent in self.sAgentSet:
                        if (random.random() <= beta1_B):
                            newI_withB_list.add(agent)
                    if (agent in self.i_A_AgentSet):
                        if(random.random() <= beta1_B*epsilon_A):
                            newI_withB_list.add(agent)
            #STANDARD CONTAGION from iAgent in AB
            #we need to loop over the agents who are currently in the state AB
            for iAgent in self.i_AB_AgentSet:
                #expose their network neighbors
                for agent in self.neighbors_dict[iAgent]:
                    if agent in self.sAgentSet:
                        if (random.random() <= beta1_B):
                            newI_withB_list.add(agent)
                    if (agent in self.i_A_AgentSet):
                        if(random.random() <= beta1_B*epsilon_A):
                            newI_withB_list.add(agent)              
            #TRIANGLE CONTAGION
            for triangle in self.triangles_list:
                n1, n2, n3 = triangle
                if (n1 in self.i_B_AgentSet) or (n1 in self.i_AB_AgentSet):
                    if (n2 in self.i_B_AgentSet) or (n2 in self.i_AB_AgentSet):
                        if n3 in self.sAgentSet: #n1, n2 infected (B or AB) n3 suscectible
                            if (random.random() <= beta2_B):
                                newI_withB_list.add(n3)
                        if (n3 in self.i_A_AgentSet):
                            if(random.random()<= beta2_B*epsilon_A):
                                newI_withB_list.add(agent) #n1, n2 infected n3 in state A
                    else: #n2 in state A or suscectible
                        if ((n3 in self.i_B_AgentSet) or (n3 in self.i_AB_AgentSet)) and (n2 in self.sAgentSet):
                            if (random.random() <= beta2_B): #n1, n3 infected (B or AB) n2 suscectible
                                newI_withB_list.add(n2)
                        if ((n3 in self.i_B_AgentSet) or (n3 in self.i_AB_AgentSet)) and (n2 in self.i_A_AgentSet):
                            if (random.random() <= beta2_B*epsilon_A): #n1, n3 infected (B or AB) n2 in state A
                                newI_withB_list.add(n2)
                else: #n1 in state A or suscectible
                    if ((n2 in self.i_B_AgentSet) or (n2 in self.i_AB_AgentSet)) and ((n3 in self.i_B_AgentSet) or (n3 in self.i_AB_AgentSet)) and (n1 in self.sAgentSet):
                        #n2, n3 infected (B or AB) and n1 suscectible
                        if (random.random() <= beta2_B):
                            newI_withB_list.add(n1)
                    if ((n2 in self.i_B_AgentSet) or (n2 in self.i_AB_AgentSet)) and ((n3 in self.i_B_AgentSet) or (n3 in self.i_AB_AgentSet)) and (n1 in self.i_A_AgentSet):
                        #n2, n3 infected (B or AB) and n1 in state A
                        if (random.random() <= beta2_B*epsilon_A):
                            newI_withB_list.add(n1)
            
            #Update only now the nodes that have been infected
            for n_to_infect in newI_withA_list:
                self.infectAgent_withA(n_to_infect)
            
            for n_to_infect in newI_withB_list:
                self.infectAgent_withB(n_to_infect)
            
            
            #for recoveries from A
            newR_fromA_list = set()            
            for recoverAgent in self.i_A_AgentSet:
                    #if the agent has just been infected it will not recover this time
                if recoverAgent in newI_withA_list:
                    continue
                else:
                    if (random.random() <= mu_A):
                        newR_fromA_list.add(recoverAgent)
            for recoverAgent in self.i_AB_AgentSet:
                if recoverAgent in newI_withA_list:
                    continue
                else:
                    if (random.random() <= mu_A):
                        newR_fromA_list.add(recoverAgent)                

            #for recoveries from B
            newR_fromB_list = set()
            for recoverAgent in self.i_B_AgentSet:
                    #if the agent has just been infected it will not recover this time
                if recoverAgent in newI_withB_list:
                    continue
                else:
                    if (random.random() <= mu_B):
                        newR_fromB_list.add(recoverAgent)
            for recoverAgent in self.i_AB_AgentSet:
                if recoverAgent in newI_withB_list:
                    continue
                else:
                    if (random.random() <= mu_B):
                        newR_fromB_list.add(recoverAgent)

            #Update only now the nodes that have recovered
            for n_to_recover in newR_fromA_list:
                self.recoverAgent_fromA(n_to_recover)            
            for n_to_recover in newR_fromB_list:
                self.recoverAgent_fromB(n_to_recover)
            
            #then track the number of individuals in each state
            self.i_A_List.append(len(self.i_A_AgentSet))
            self.i_B_List.append(len(self.i_B_AgentSet))
            self.i_AB_List.append(len(self.i_AB_AgentSet))
            self.S_List.append(len(self.sAgentSet))
            
            #increment the time
            self.t += 1

        return self.i_A_List, self.i_B_List, self.i_AB_List, self.S_List
    
    def get_stationary_rho_A(self, normed=True, last_k_values = 100):
        i = self.i_A_List
        if len(i)==0:
            return 0
        if normed:
            i = 1.*np.array(i)/self.N
        if i[-1]==1:
            return 1
        elif i[-1]==0:
            return 0
        else:
            avg_i = np.mean(i[-last_k_values:])
            avg_i = np.nan_to_num(avg_i) #if there are no infected left nan->0
        return avg_i

    def get_stationary_rho_B(self, normed=True, last_k_values = 100):
        i = self.i_B_List
        if len(i)==0:
            return 0
        if normed:
            i = 1.*np.array(i)/self.N
        if i[-1]==1:
            return 1
        elif i[-1]==0:
            return 0
        else:
            avg_i = np.mean(i[-last_k_values:])
            avg_i = np.nan_to_num(avg_i) #if there are no infected left nan->0
            return avg_i
    
    def get_stationary_rho_AB(self, normed=True, last_k_values = 100):
        i = self.i_AB_List
        if len(i)==0:
            return 0
        if normed:
            i = 1.*np.array(i)/self.N
        if i[-1]==1:
            return 1
        elif i[-1]==0:
            return 0
        else:
            avg_i = np.mean(i[-last_k_values:])
            avg_i = np.nan_to_num(avg_i) #if there are no infected left nan->0
            return avg_i
    

def generate_my_simplicial_complex_d2(N,p1,p2):
    #I first generate a standard ER graph with edges connected with probability p1
    G = nx.fast_gnp_random_graph(N, p1, seed=None)
    if not nx.is_connected(G):
        giant = max(nx.connected_components(G), key=len)
        G = G.subgraph(giant).copy()
        print ('not connected, but GC has order ', G.order(), 'and size', G.size())
    triangles_list = []
    
    #Now I run over all the possible combinations of three elements:
    for tri in combinations(list(G.nodes()),3):
        #And I create the triangle with probability p2
        if random.random() <= p2:
            #I close the triangle.
            triangles_list.append(tri)
            
            #Now I also need to add the new links to the graph created by the triangle
            G.add_edge(tri[0], tri[1])
            G.add_edge(tri[1], tri[2])
            G.add_edge(tri[0], tri[2])
             
    #Creating a dictionary of neighbors
    node_neighbors_dict = {}
    for n in G.nodes():
        node_neighbors_dict[n] = G[n].keys()
                
    print (len(triangles_list), 'triangles created. Size now is', G.size())

    return node_neighbors_dict, triangles_list

def get_p1_and_p2(k1,k2,N):
    p2 = (2.*k2)/((N-1.)*(N-2.))
    p1 = (k1 - 2.*k2)/((N-1.)- 2.*k2)
    if (p1>=0) and (p2>=0):
        return p1, p2
    else:
        raise ValueError('Negative probability!')
        


        
def run_one_simulation(args):
    
    it_num, N, p1, p2, lambda1s_A, lambda1s_B, lambdaD_A_target, lambdaD_B_target, I_A_percentage, I_B_percentage, I_AB_percentage, epsilon_A, epsilon_B, t_max, mu_A, mu_B  = args 
    
    print('simulation', it_num, 'has started')
    node_neighbors_dict, triangles_list = generate_my_simplicial_complex_d2(N,p1,p2)

    real_k = 1.*sum([len(v) for v  in node_neighbors_dict.values()])/len(node_neighbors_dict)
    real_kD = 3.*len(triangles_list)/len(node_neighbors_dict)
    
    beta1s_A = []
    for lambda1_A in lambda1s_A:
        beta1_A = 1.*(mu_A/real_k)*lambda1_A
        beta1s_A.append(beta1_A)

    beta2_A = 1.*(mu_A/real_kD)*lambdaD_A_target

    beta1s_B = []
    for lambda1_B in lambda1s_B:
        beta1_B = 1.*(mu_B/real_k)*lambda1_B
        beta1s_B.append(beta1_B)

    beta2_B = 1.*(mu_B/real_kD)*lambdaD_B_target    
    rhos_A = []
    rhos_B = []
    rhos_AB = []
    
    for (beta1_A, beta1_B) in zip(beta1s_A, beta1s_B):
        mySimplagionModel = SimplagionModel(node_neighbors_dict, triangles_list, I_A_percentage, I_B_percentage, I_AB_percentage)
        mySimplagionModel.initial_setup(fixed_nodes_to_infect_withA = mySimplagionModel.infected_this_setup_withA, fixed_nodes_to_infect_withB=mySimplagionModel.infected_this_setup_withB, fixed_nodes_to_infect_withAB=mySimplagionModel.infected_this_setup_withAB);
        results = mySimplagionModel.run(t_max, beta1_A, beta2_A, mu_A, beta1_B, beta2_B, mu_B, epsilon_A, epsilon_B)
        rho_A = mySimplagionModel.get_stationary_rho_A(normed=True, last_k_values = 100)
        rhos_A.append(rho_A)
        rho_B = mySimplagionModel.get_stationary_rho_B(normed=True, last_k_values = 100)
        rhos_B.append(rho_B)
        rho_AB = mySimplagionModel.get_stationary_rho_AB(normed=True, last_k_values = 100)
        rhos_AB.append(rho_AB)

    return rhos_A, rhos_B, rhos_AB, real_k, real_kD

def parse_results(results):
    
    rhos_A_array, rhos_B_array, rhos_AB_array, real_k_list, real_kD_list = [], [], [], [], []
    
    for rhos_A, rhos_B, rhos_AB, real_k, real_kD in results:
        real_k_list.append(real_k)
        real_kD_list.append(real_kD)
        rhos_A_array.append(rhos_A)
        rhos_B_array.append(rhos_B)
        rhos_AB_array.append(rhos_AB)
        
        
    rhos_A_array = np.array(rhos_A_array)
    rhos_B_array = np.array(rhos_B_array)
    rhos_AB_array = np.array(rhos_AB_array)
    real_kD_list = np.array(real_kD_list)
    real_k_list = np.array(real_k_list)
    
    avg_kD = real_kD_list.mean(axis=0)
    avg_k = real_k_list.mean(axis=0)

    avg_rhos_A = np.mean(rhos_A_array, axis=0)
    avg_rhos_B = np.mean(rhos_B_array, axis=0)
    avg_rhos_AB = np.mean(rhos_AB_array, axis=0)
        
    return avg_rhos_A, avg_rhos_B, avg_rhos_AB, avg_k, avg_kD
        