import os
import yaml
from CPAC.utils import Configuration
from CPAC.pipeline.cpac_pipeline import prep_workflow
from CPAC.pipeline.cpac_runner import build_strategies
import unittest
import networkx as nx

        
class test_graph_of_pipeline():
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

    def __init__(self):
        
        # load stored graph
        
       
        # get the subject list
        subListPath = "/home/rtungaraza/CPAC_testData/CPAC_subjectList.yml"
        self.sublist = yaml.load(open(os.path.realpath(subListPath), 'r'))    
            
        
        # get pipeline configuration file
        configFile = "/home/rtungaraza/CPAC_testData/"+\
                        "pipeline_config_test_8subjects_pipeline.yml"
        self.c = Configuration(yaml.load(open(os.path.realpath(configFile),\
                                'r')))

        # build strategies 
        self.strategies = None          
        if (1 in self.c.runNuisance) or (self.c.Corrections != None):
            self.strategies = sorted(build_strategies(self.c))
        
        # Run the pipeline building 
        self.workflow = prep_workflow(self.sublist[0], self.c, \
                                        self.strategies, 0)
    
    def concat_edges(self,edge):
        print edge
        return "<-->".join([edge[0].name,edge[1].name])    
    
    def test_duplicate_edges(self):
        items = self.workflow._graph.edges()
        items = map(list, items)
        print map(self.concat_edges, items)
        
    def test_pipeline_graph(self):
        for n,nbrs in self.workflow._graph.adjacency_iter():
            for nbr,eattr in nbrs.items():
                data=eattr['connect']
                print("%s, %s, %s\n" % (n,nbr,data))

##    def test_pipeline_number_of_nodes(self):
##        print "Number of nodes: ",self.workflow._graph.number_of_nodes()
##        
##    def test_pipeline_number_of_edges(self):
##        pass
##        
##    def test_pipeline_edges(self):
##        pass
##        
##    def test_pipeline_list_of_nodes(self):
##        print "List of nodes: ",self.workflow._graph.nodes()
        
    







me = test_graph_of_pipeline()
me.test_duplicate_edges()

























        
