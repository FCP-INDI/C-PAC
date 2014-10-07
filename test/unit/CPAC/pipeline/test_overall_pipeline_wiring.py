import os
import yaml
from CPAC.utils import Configuration
from CPAC.pipeline.cpac_pipeline import prep_workflow
from CPAC.pipeline.cpac_runner import build_strategies
from operator import itemgetter
import unittest
import pickle
import matplotlib.pyplot as plt
import networkx as nx

        
class TestPipelineGraph(unittest.TestCase):
    """
    Objective: This class tests whether any changes to cpac_pipeline.py
               have globally changed the expected graph resulting from the 
               pre-determined pipeline configuration file and subject list
    Author:    Rosalia Tungaraza
    Date:      September 30, 2014

    Input:     subListPath -> a predetermined subject list
               configFile  -> a predetermined pipeline configuration file 
            
    Output:    true/false depending on whether the new created nipype 
               graph/pipeline is exactly the same as the expected graph
    """

    def setUp(self):
        
        # load the gold standard graph  
        edgelistFile = "pipeline_edge_list.p"
        self.edge_list=pickle.load(open(edgelistFile, "rb" ))         

        edgePropFile = "pipeline_edgeProperties_list.p"
        self.edgeProperties = pickle.load(open(edgePropFile,"rb"))           

        graphOverviewFile = "pipeline_graph_overview.p"
        self.g_overview=pickle.load(open(graphOverviewFile,"rb"))           

        subFile = "subject_list_4_unittest.yml"
        self.sublist = yaml.load(open(subFile, 'r'))  

        configFile = "pipeline_4_unittest.yml"
        self.c = Configuration(yaml.load(open(configFile,'r')))

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
                                            
    def test_edge_list(self):
        new_edge_list = self.workflow._graph.edges()
        new_edge_list = map(self.extract_name, new_edge_list)
        new_edge_list = sorted(new_edge_list, key=itemgetter(0,1))
        self.assertSequenceEqual(new_edge_list,self.edge_list,\
                    msg="Edge list differ between the two graphs")
        
    def test_edge_properties(self):        
        new_edgeProperties=[]       
        for n,nbrs in self.workflow._graph.adjacency_iter():
            for nbr,eattr in nbrs.items():
                data=eattr['connect']
                new_edgeProperties.append([n.name,nbr.name,repr(data)])        
        new_edgeProperties = sorted(new_edgeProperties, key=itemgetter(0,1,2))        
        self.assertSequenceEqual(new_edgeProperties,self.edgeProperties,\
                    msg="Edge properties differ between the two graphs")
        
    def test_list_of_nodes(self): 
        nodes=sorted(map(self.extract_name,self.workflow._graph.nodes()))        
        self.assertSequenceEqual(nodes, self.g_overview["list_of_nodes"], \
                        msg="List of nodes differ between the two graphs")
            
    def test_total_nodes(self):
        self.assertEqual(self.g_overview["total_nodes"],\
                    self.workflow._graph.number_of_nodes(),\
                    msg="Total Number of nodes differ between the two graphs")
        
    def test_total_edges(self):
        self.assertEqual(self.g_overview["total_edges"],\
                    self.workflow._graph.number_of_edges(),\
                    msg="Total Number of edges differ between the two graphs")
            
    # def test_graph_draw(self):
        

#------------------------------------------------------------------------------
# run test
if __name__ == '__main__':
    unittest.main()


















        
