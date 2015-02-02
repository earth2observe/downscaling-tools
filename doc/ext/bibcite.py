from docutils import nodes
from docutils.parsers.rst.directives import unchanged

from sphinx.util.compat import Directive

def setup(app):
    app.add_node(cite,
                 html=(cite.visit_node, cite.depart_node),
                 latex=(cite.visit_node, cite.depart_node),
                 text=(cite.visit_node, cite.depart_node))
    app.add_generic_role('cite', cite)
    app.add_node(citep,
                 html=(citep.visit_node, citep.depart_node),
                 latex=(citep.visit_node, citep.depart_node),
                 text=(citep.visit_node, citep.depart_node))
    app.add_generic_role('citep', citep)
    app.add_node(bibliography,
                 html=(bibliography.visit_node, bibliography.depart_node),
                 latex=(bibliography.visit_node, bibliography.depart_node),
                 text=(bibliography.visit_node, bibliography.depart_node))
    app.add_directive('bibliography', bibliographyDirective)

class cite(nodes.Inline, nodes.TextElement):

    @staticmethod
    def visit_node(translator, node):
        if node.children:
            if translator.builder.name == 'latex':
                translator.body.append(u"\\%s{" % node.tagname)
            elif translator.builder.name == 'html':
                pre = u"<em>(" if node.tagname == 'citep' else u'<span class="citation">'
                translator.body.append(pre)

    @staticmethod
    def depart_node(translator, node):
        if node.children:
            if translator.builder.name == 'latex':
                translator.body.append(u"} ")
            elif translator.builder.name == 'html':
                post = u")</em>" if node.tagname == 'citep' else u"</span>"
                translator.body.append(post)

class citep(cite):
    pass

class bibliography(nodes.paragraph):
    """TODO: currently not used. Remove and keep ?"""

    @staticmethod
    def visit_node(translator, node):
        pass

    @staticmethod
    def depart_node(translator, node):
        pass

class bibliographyDirective(Directive):

    has_content = False
    required_arguments = 1
    option_spec = {'style': unchanged}

    def run(self):
        style = self.options['style']
        bibdb = self.arguments[0]
        bibfile = bibdb + '.bib'
        self.state.document.settings.env.config.latex_additional_files.append(bibfile)
        attributes = {'format': 'latex'}
        style_node = nodes.raw('', "\\bibliographystyle{%s}" % style,
                               **attributes)
        bibdb_node = nodes.raw('', "\\bibliography{%s}" % bibdb,
                               **attributes)
        return [style_node, bibdb_node]

