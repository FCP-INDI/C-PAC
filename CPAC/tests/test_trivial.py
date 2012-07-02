from nose.tools import ok_, eq_

def test_b():
    """
    Raw, unparented test.
    """
    assert 'b' == 'b'

def test_1_and_1():
    assert 1+1 == 2
    
def test_sum():
    eq_(1+1,2)

def test_failing_compare():
    ok_(2>3, 'Expected failure')
    
