from ddt import ddt, data, unpack
from django.test import TestCase
from commons.models import  User


@ddt
class TestUser(TestCase):

    def test_CRUD(self):
        email = 'test@example.com'
        # create
        res = User.objects.create(email=email, password="test")
        assert str(res) == email

        # get
        res = User.objects.get(email=email)
        assert str(res) == email

        #update
        res = User.objects.filter(email=email).update(password="test2")
        assert res == 1
        res = User.objects.filter(email='wrong').update(password="test3")
        assert res == 0

        #delete
        res = User.objects.filter(email=email).delete()
        assert res[0] == 1
    
    @data(
        ['test@example.com', 'test', None, False],
        ['test@example.com', 'test', True, True],
    )
    @unpack
    def test_create_user(self, email, password, is_staff, expect):
        res = User.objects.create_user(email, password, is_staff)
        assert res.is_staff == expect

    @data(
        ['test@example.com', 'test', True],
    )
    @unpack
    def test_create_staff(self, email, password, expect):
        res = User.objects.create_staff(email, password)
        assert res.is_staff == expect
        print(f"##{res.is_staff}##")   
    