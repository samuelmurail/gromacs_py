
# Remove all script from 'gromacs_py' directory:
# Remove all script from 'notebook' directory:
# Remove all script from 'ci-helpers' directory: (Travis issue)

def pytest_ignore_collect(path):
    path_list = str(path).split('/')
    if path_list[-1]!='__init__.py' and path_list[-1][-3:]=='.py' and path_list[-2] in ['gromacs_py', 'notebook', 'ci-helpers']:
        print('Ignore :',str(path))
        return True
