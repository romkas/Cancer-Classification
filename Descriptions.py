#covariates2002 = [
#    ('Rand1', 'null'),
#    ('Tod', 'null'),
#    ('Time', 'null'),
#    ('Sex', 'cat'),
#    ('Immun', 'cat'),
#    ('CNS', 'cat'),
#    ('Mediastinum', 'cat'),
#    ('Age', 'num'),
#    ('Leuc', 'num'),
#    ('Leber', 'num'),
#    ('Milz', 'num'),
#]
covariates2002 = ['Rand',
                  'Tod',
                  'Time',
                  'Sex',
                  'Immun',
                  'CNS',
                  'Mediastinum',
                  'Age',
                  'Leuc',
                  'Leber',
                  'Milz']
min_sub_size = 100
min_time_research = 50
num_split_at_level = 5

"""
covLevels_2002 = [
    [],
    [],
    [],
    [1, 2],
    [1, 2, 3, 5, 6, 7],
    [1, 2],
    [1, 2],
    [3, 7, 12, 16],
    [10, 15, 20, 25, 30, 35, 40],
    [2, 4, 6, 8, 10, 12, 14],
    [2, 4, 6, 8, 10, 12, 14, 16, 18],
]

covariates_2008 = [
    ('Rand1', 'null'),
    ('Tod', 'null'),
    ('Time', 'null'),
    ('Sex', 'cat'),
    ('Immun1-3', 'cat'),
    ('Immun10', 'cat'),
    ('Immun12', 'cat'),
    ('Immun13', 'cat'),
    ('Immun11', 'cat'),
    ('Immun14-20', 'cat'),
    ('CNS', 'cat'),
    ('Mediastinum', 'cat'),
    ('Age', 'num'),
    ('Leuc', 'num'),
    ('Leber', 'num'),
    ('Milz', 'num'),
]

covLevels_2008 = [
    [],
    [],
    [],
    [1, 2],
    [0, 1],
    [0, 1],
    [0, 1],
    [0, 1],
    [0, 1],
    [0, 1],
    [1, 2, 3],
    [1, 2],
    [3, 7, 12, 16],
    [10, 15, 20, 25, 30, 35, 40],
    [2, 4, 6, 8, 10, 12, 14],
    [2, 4, 6, 8, 10, 12, 14, 16, 18],
]

minSubSize = 100
treeHeight_2002 = len(covariates_2002)
treeHeight_2008 = len(covariates_2008)
numBestCov = 5
minTimeResearch = 50

max_len_cov_name = len('Mediastinum')
"""