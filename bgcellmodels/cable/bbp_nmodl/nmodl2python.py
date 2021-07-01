"""
Converting NMODL to Python code using example copied from 
https://bluebrain.github.io/nmodl/html/notebooks/nmodl-python-tutorial.html
"""

import nmodl.dsl as nmodl
from nmodl.dsl import ast
from nmodl.dsl import visitor

class PyGenerator(visitor.AstVisitor):
    def __init__(self):
        visitor.AstVisitor.__init__(self)
        self.pycode = ''
        self.indent = 0
        self.func_name = ""

    def visit_function_block(self, node):
        params = []
        self.func_name = node.get_node_name()
        for p in node.parameters:
            params.append(p.get_node_name())
        params_str = ", ".join(params)
        self.pycode += f"def {node.get_node_name()}({params_str}):\n"
        node.visit_children(self)

    def visit_statement_block(self, node):
        self.indent += 1
        node.visit_children(self)
        self.indent -= 1

    def visit_expression_statement(self, node):
        self.pycode += " "*4*self.indent
        expr = node.expression
        if type(expr) is ast.BinaryExpression and expr.op.eval() == "=":
            rhs = expr.rhs
            lhsn = expr.lhs.name.get_node_name()
            if lhsn == self.func_name:
                self.pycode += "return "
                rhs.accept(self)
            else:
                node.visit_children(self)
        else:
            node.visit_children(self)
        self.pycode += "\n"


    def visit_if_statement(self, node):
        self.pycode += " "*4*self.indent + "if "
        node.condition.accept(self)
        self.pycode += ":\n"
        node.get_statement_block().accept(self)
        for n in node.elseifs:
            n.accept(self)
        if node.elses:
            node.elses.accept(self)

    def visit_else_statement(self, node):
        self.pycode += " "*4*self.indent + "else:\n"
        node.get_statement_block().accept(self)


    def visit_binary_expression(self, node):
        lhs = node.lhs
        rhs = node.rhs
        op = node.op.eval()
        if op == "^":
            self.pycode += "pow("
            lhs.accept(self)
            self.pycode += ", "
            rhs.accept(self)
            self.pycode += ")"
        else:
            lhs.accept(self)
            self.pycode += f" {op} "
            rhs.accept(self)

    def visit_var_name(self, node):
        self.pycode += node.name.get_node_name()

    def visit_integer(self, node):
        self.pycode += nmodl.to_nmodl(node)

    def visit_double(self, node):
        self.pycode += nmodl.to_nmodl(node)

if __name__ == '__main__':
    # Parse NMODL code
    driver = nmodl.NmodlDriver()
    mfunc_ast = driver.parse_string(mfunc_src)

    # Convert to Python code
    pygen = PyGenerator()
    pygen.visit_program(mfunc_ast)
    print(pygen.pycode)