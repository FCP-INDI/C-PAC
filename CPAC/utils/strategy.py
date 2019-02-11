import os
import six
import warnings
import logging

logger = logging.getLogger('workflow')


class Strategy(object):

    def __init__(self):
        self.resource_pool = {}
        self.leaf_node = None
        self.leaf_out_file = None
        self.name = []

    def append_name(self, name):
        self.name.append(name)

    def get_name(self):
        return self.name

    def set_leaf_properties(self, node, out_file):
        self.leaf_node = node
        self.leaf_out_file = out_file

    def get_leaf_properties(self):
        return self.leaf_node, self.leaf_out_file

    def get_resource_pool(self):
        return self.resource_pool

    def get_nodes_names(self):
        pieces = [n.split('_') for n in self.name]
        assert all(p[-1].isdigit() for p in pieces)
        return ['_'.join(p[:-1]) for p in pieces]

    def get_node_from_resource_pool(self, resource_key):
        try:
            return self.resource_pool[resource_key]
        except:
            logger.error('No node for output: %s', resource_key)
            raise

    def update_resource_pool(self, resources, override=False):
        for key, value in resources.items():
            if key in self.resource_pool and not override:
                raise Exception(
                    'Key %s already exists in resource pool, '
                    'replacing with %s ' % (key, value)
                )
            self.resource_pool[key] = value

    def __getitem__(self, resource_key):
        assert isinstance(resource_key, six.string_types)
        try:
            return self.resource_pool[resource_key]
        except:
            logger.error('No node for output: %s', resource_key)
            raise

    def __contains__(self, resource_key):
        assert isinstance(resource_key, six.string_types)
        return resource_key in self.resource_pool

    def fork(self):
        fork = Strategy()
        fork.resource_pool = dict(self.resource_pool)
        fork.leaf_node = self.leaf_node
        fork.out_file = str(self.leaf_out_file)
        fork.leaf_out_file = str(self.leaf_out_file)
        fork.name = list(self.name)
        return fork
