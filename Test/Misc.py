'''
Misc test.
'''

__all__=['test_misc']

def test_misc(arg):
    if arg in ('tree','misc','all'):
        from HamiltonianPy.Misc.test import test_tree
        test_tree()
    if arg in ('linalg','misc','all'):
        from HamiltonianPy.Misc.test import test_linalg
        test_linalg()
