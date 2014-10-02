import os
import yaml
from CPAC.utils import Configuration
from CPAC.pipeline.cpac_pipeline import prep_workflow
from CPAC.pipeline.cpac_runner import build_strategies
from operator import itemgetter
import matplotlib.pyplot as plt
import networkx as nx
import unittest
import pickle


        
class GoldStandardGraph(object):
    """
    Objective:  This class tests creates the graph that will be stored in the
                test/unit/ directory and will be used as gold standard against
                all graphs that will be generated when modifying the
                cpac_pipeline.py code. 
    Author:    Rosalia Tungaraza
    Date:      September 30, 2014

    Input:     subFile -> a predetermined subject list
               configFile  -> a predetermined pipeline configuration file 
            
    Output:    true/false depending on whether the new created nipype 
               graph/pipeline is exactly the same as the expected graph
    """

    def __init__(self):
        
        # get the subject list
        subFile = "subject_list_4_unittest.yml"
        self.sublist = yaml.load(open(os.path.realpath(subFile), 'r'))    
            
        
        # get pipeline configuration file
        configFile = "pipeline_4_unittest.yml"
        self.c = Configuration(yaml.load(open(os.path.realpath(configFile),\
                                'r')))

        # build strategies 
        self.strategies = None          
        if (1 in self.c.runNuisance) or (self.c.Corrections != None):
            self.strategies = sorted(build_strategies(self.c))
        
        # Run the pipeline building 
        self.workflow = prep_workflow(self.sublist[0], self.c, \
                                        self.strategies, 0)
    def extract_name(self,node):
        if isinstance(node, (list, tuple)):
            return [node[0].name,node[1].name]            
        else:
            return node.name
        
    def generate_edge_list(self):
        self.edges = self.workflow._graph.edges()
        self.edges = map(self.extract_name, self.edges)
        self.edges = sorted(self.edges, key=itemgetter(0,1))
        try: 
            pickle.dump( self.edges, open( "pipeline_edge_list.p", "wb" ) )           
        except:
            print "Couldn't pickle edge list ...\n"
        
    def generate_edge_properties(self): 
        self.edgeProperties=[]       
        for n,nbrs in self.workflow._graph.adjacency_iter():
            for nbr,eattr in nbrs.items():
                data=eattr['connect']
                self.edgeProperties.append([n.name,nbr.name,repr(data)])
        
        self.edgeProperties=sorted(self.edgeProperties, key=itemgetter(0,1,2))
        try: 
            pickle.dump(self.edgeProperties, \
                        open( "pipeline_edgeProperties_list.p", "wb" ) )           
        except:
            print "Couldn't pickle edgeProperties list ...\n"      
        
    def generate_graph_overview(self):
        self.g_overview=dict()
        self.g_overview["list_of_nodes"] = \
                    sorted(map(self.extract_name,self.workflow._graph.nodes()))
        self.g_overview["total_nodes"] = self.workflow._graph.number_of_nodes()
        self.g_overview["total_edges"] = self.workflow._graph.number_of_edges()
        try: 
            pickle.dump(self.g_overview, \
                        open( "pipeline_graph_overview.p", "wb" ) )           
        except:
            print "Couldn't pickle graph overview dictionary ...\n" 
    
    def find_root_node(self,graph):  
        roots=[node for node,degree in graph.in_degree().items() if degree==0]
        return sorted(roots)
                        
          
    def write_graph(self):
        tgraph = nx.DiGraph()
        tgraph.add_nodes_from(self.g_overview["list_of_nodes"])
        tgraph.add_edges_from(self.edges)             
        nx.write_graphml(tgraph, "test.graphml")  
        # run "dot -Tpng test.dot >test.png" to create png image
        nx.write_dot(tgraph,'test.dot')     
    
#------------------------------------------------------------------------------
# Testing the code
if __name__ == "__main__":
    g = GoldStandardGraph()
    g.generate_edge_list()
    g.generate_edge_properties()
    g.generate_graph_overview()
    g.write_graph()






































