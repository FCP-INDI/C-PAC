class SchemaError(Exception):
    """Exception raised for errors in the schema."""

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class NodeBlockError(Exception):
    """Exception raised for errors in the node block."""

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)
