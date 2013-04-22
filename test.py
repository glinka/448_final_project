def test(v):
    reload(v)
    A, Opns, p = v.vote(n=100)
    v.checkDegrees(A, p)
