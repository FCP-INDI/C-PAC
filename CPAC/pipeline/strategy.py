


class strategy:

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

    def get_node_from_resource_pool(self, resource_key):
        try:
            if resource_key in self.resource_pool:
                return self.resource_pool[resource_key]
        except:
            logger.info('no node for output: ')
            logger.info(resource_key)
            raise

    def update_resource_pool(self, resources):
        for key, value in resources.items():
            if key in self.resource_pool:
                logger.info('Warning key %s already exists in resource' \
                        ' pool, replacing with %s ' % (key, value))

            self.resource_pool[key] = value