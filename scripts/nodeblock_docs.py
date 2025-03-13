import ast
import json
import os
from pathlib import Path

class NodeBlockVisitor(ast.NodeVisitor):
    def __init__(self):
        self.nodeblocks = []
        self.current_file = None

    def safe_eval(self, node):
        """Safely evaluate AST nodes, handling special cases."""
        if isinstance(node, ast.Constant):
            return node.value
        elif isinstance(node, ast.List):
            return [self.safe_eval(elt) for elt in node.elts]
        elif isinstance(node, ast.Tuple):
            return tuple(self.safe_eval(elt) for elt in node.elts)
        elif isinstance(node, ast.Starred):
            # For starred expressions, we'll return a placeholder
            return f"*{self.get_starred_expr(node.value)}"
        elif isinstance(node, ast.Name):
            return f"@{node.id}"  # Return variable names with @ prefix to distinguish them
        elif isinstance(node, ast.Attribute):
            # Handle attribute access (e.g., module.attribute)
            return f"@{self.get_attribute_chain(node)}"
        else:
            return f"!UNSUPPORTED_TYPE:{type(node).__name__}!"

    def get_starred_expr(self, node):
        """Get a string representation of what's being starred."""
        if isinstance(node, ast.Name):
            return node.id
        elif isinstance(node, ast.Attribute):
            return self.get_attribute_chain(node)
        return f"unknown_starred_{type(node).__name__}"

    def get_attribute_chain(self, node):
        """Get the full chain of attributes (e.g., 'module.attribute')."""
        parts = []
        current = node
        while isinstance(current, ast.Attribute):
            parts.append(current.attr)
            current = current.value
        if isinstance(current, ast.Name):
            parts.append(current.id)
        return '.'.join(reversed(parts))

    def visit_FunctionDef(self, node):
        # Check if the function has a decorator that matches @nodeblock
        for decorator in node.decorator_list:
            if isinstance(decorator, ast.Call):
                if isinstance(decorator.func, ast.Name) and decorator.func.id == 'nodeblock':
                    # Extract decorator arguments
                    args = {}
                    for kw in decorator.keywords:
                        # Use safe_eval instead of ast.literal_eval
                        args[kw.arg] = self.safe_eval(kw.value)
                    
                    # Get function docstring
                    docstring = ast.get_docstring(node)
                    
                    # Get source code lines
                    source_lines = []
                    for i in range(node.lineno - 1, node.end_lineno):
                        source_lines.append(self.current_source[i])
                    source_code = '\n'.join(source_lines)

                    # Create nodeblock info
                    nodeblock_info = {
                        'name': node.name,
                        'file': str(self.current_file),
                        'line_number': node.lineno,
                        'decorator_args': args,
                        'docstring': docstring,
                        'source_code': source_code,
                    }
                    self.nodeblocks.append(nodeblock_info)

def find_nodeblocks(root_dir):
    visitor = NodeBlockVisitor()
    root_path = Path(root_dir)

    # Walk through all Python files
    for python_file in root_path.rglob('*.py'):
        try:
            with open(python_file, 'r', encoding='utf-8') as f:
                source = f.read()
                visitor.current_source = source.splitlines()
                visitor.current_file = python_file.relative_to(root_path)
                tree = ast.parse(source)
                visitor.visit(tree)
        except Exception as e:
            print(f"Error processing {python_file}: {e}")

    return visitor.nodeblocks

def main():
    # Assuming you're running this from the root of the C-PAC repository
    root_dir = '.'  # or specify the full path to C-PAC repository
    nodeblocks = find_nodeblocks(root_dir)

    # Save to JSON file
    output_file = 'nodeblock_index.json'
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(nodeblocks, f, indent=2, ensure_ascii=False)

    print(f"Found {len(nodeblocks)} nodeblocks")
    print(f"Index saved to {output_file}")

if __name__ == '__main__':
    main()
