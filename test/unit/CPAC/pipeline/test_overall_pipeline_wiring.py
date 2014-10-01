import os
import yaml
from CPAC.utils import Configuration
from CPAC.pipeline.cpac_pipeline import prep_workflow
from CPAC.pipeline.cpac_runner import build_strategies
import unittest
import networkx as nx
import pickle

        
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
        try: 
            self.edge_list=pickle.load(open("pipeline_edge_list.p", "rb" ))           
        except:
            print "Couldn't open pickled edge list ...\n"
            
        try: 
            self.edgeProperties = \
                pickle.load(open("pipeline_edgeProperties_list.p","rb"))           
        except:
            print "Couldn't open pickled edgeProperties list ...\n"  
        
        try: 
            self.g_overview=pickle.load(open("pipeline_graph_overview.p","rb"))           
        except:
            print "Couldn't open pickled graph overview dictionary ...\n"   
              
        # get the subject list
        subFile = "subject_list_4_unittest.yml"
        try:
            self.sublist = yaml.load(open(subFile), 'r')  
        except:
            print "Couldn't open file with SUBJECT LIST ...\n"              
        
        # get pipeline configuration file
        configFile = "pipeline_4_unittest.yml"
        try:
            self.c = Configuration(yaml.load(open(configFile),'r'))
        except:
            print  "Couldn't open file with PIPELINE CONFIGURATION ...\n" 

        # build strategies 
        self.strategies = None          
        if (1 in self.c.runNuisance) or (self.c.Corrections != None):
            self.strategies = sorted(build_strategies(self.c))
        
        # Run the pipeline building 
        self.workflow = prep_workflow(self.sublist[0], self.c, \
                                        self.strategies, 0)
                                        
    def test_edge_list(self):
        new_edge_list = self.workflow._graph.edges()
        new_edge_list = sorted(new_edge_list, key=itemgetter(0))
        self.assertEqual(new_edge_list,self.edge_list,\
                    msg="Edge list differ between the two graphs")
        
    def test_edge_properties(self):        
        new_edgeProperties=[]       
        for n,nbrs in self.workflow._graph.adjacency_iter():
            for nbr,eattr in nbrs.items():
                data=eattr['connect']
                new_edgeProperties.append([n,nbr,data])        
        new_edgeProperties = sorted(new_edgeProperties, key=itemgetter(0))        
        self.assertEqual(new_edgeProperties,self.edgeProperties,\
                    msg="Edge properties differ between the two graphs")
        
    def test_list_of_nodes(self):        
        self.assertEqual(self.g_overview["list_of_nodes"], \
                    sorted(self.workflow._graph.nodes()),\
                    msg="List of nodes differ between the two graphs")
            
    def test_total_nodes(self):
        self.assertEqual(self.g_overview["total_nodes"],\
                    self.workflow._graph.number_of_nodes(),\
                    msg="Total Number of nodes differ between the two graphs")
        
    def test_total_edges(self):
        self.assertEqual(self.g_overview["total_edges"],\
                    self.workflow._graph.number_of_edges(),\
                    msg="Total Number of edges differ between the two graphs")
            























        
