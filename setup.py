from distutils.core import setup

setup(
    name='DWSTree',
    version='0.1.0',
    author='Dylan SCwhilk',
    author_email='dylan@schwilk.org',
    packages=['dwstree'],
    scripts=['scripts/treematic.py', 
             'scripts/treeutils.py', 
             'scripts/rarefaction.py', 
             'scripts/prune2orders.py'],
    url='https://github.com/dschwilk/dwstree',
   # license='LICENSE.txt',
    description='Phylogenetic tools.',
    long_description=open('README.md').read(),
) 
