import os
import yaml
from CPAC.utils import Configuration
from CPAC.pipeline.cpac_pipeline import prep_workflow
from CPAC.pipeline.cpac_runner import build_strategies
import nose
import networkx as nx

        
def test_WholePipeline():
    """
    Objective: This function tests whether any changes to cpac_pipeline.py
               have globally changed the expected pipeline configuration 
    Author:    Rosalia Tungaraza
    Date:      September 30, 2014

    Input:     
    Output:
    """

    # get the subject list
    subListPath = "/home/rtungaraza/CPAC_testData/CPAC_subjectList.yml"
    sublist = yaml.load(open(os.path.realpath(subListPath), 'r'))    
            
        
    # get pipeline configuration file
    configFile = "/home/rtungaraza/CPAC_testData/pipeline_config_test_8subjects_pipeline.yml"
    c = Configuration(yaml.load(open(os.path.realpath(configFile), 'r')))

    # build strategies           
    if (1 in c.runNuisance) or (c.Corrections != None):
        strategies = sorted(build_strategies(c))
    else:
        strategies = None
                
    # Run the pipeline building   
           
    workflow = prep_workflow(sublist[0], c, strategies, 0)
    print workflow.base_dir
##    print "Number of nodes: ",workflow._graph.number_of_nodes()

    for n,nbrs in workflow._graph.adjacency_iter():
        for nbr,eattr in nbrs.items():
            data=eattr['connect']
            print("%s, %s, %s\n" % (n,nbr,data))
        
        
##    print "List of nodes: ",workflow._graph.nodes()
##    nx.write_adjlist(workflow._graph, "test.adjlist")
##            


test_WholePipeline()






























        
