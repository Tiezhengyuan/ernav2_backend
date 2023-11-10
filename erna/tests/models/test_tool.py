from ddt import ddt, data, unpack
import json
from django.test import TestCase, override_settings
from commons.models import Tool

@ddt
class TestTool(TestCase):
  
    def test_CRUD(self):
        # create
        res = Tool.objects.create(
            tool_name='bowtie',
            tool_path = '/mnt/tools/',
            version = '3.4.5',
        )
        assert str(res) == 'bowtie'

        # get
        res = Tool.objects.get(tool_name='bowtie')
        assert res.version == '3.4.5'

        # udpate
        res = Tool.objects.filter(tool_name='bowtie', version='3.4.5')\
            .update(version='4.1.0')
        assert res == 1
        res = Tool.objects.filter(tool_name='bowtie', version='3.4.5')\
            .update(version='4.1.0')
        assert res == 0

        # delete
        res = Tool.objects.filter(tool_name='bowtie').delete()
        assert res[0] == 1
        res = Tool.objects.filter(tool_name='bowtie').delete()
        assert res[0] == 0
    
    @data(
        # new
        [
            {
                'tool_name': 'bowtie',
                'tool_path': 'bowtie/4.1/bowtie',
                'version':'4.1',
                'default_parameters': json.dumps({'a':4, 'b':"abc"})
            },
            'bowtie',
            json.dumps({"a": 4, "b": "abc"}),
        ],
    )
    @unpack
    def test_add_tool(self, input, exepct_name, expect_params):
        res = Tool().add_tool(input)
        assert res.tool_name == exepct_name
        assert res.default_parameters == expect_params
        
    
    def test_tools_list(self):
        res = Tool.objects.tools_list()
        assert res == {}

        Tool.objects.create(
            tool_name='Bowtie',
            tool_path = 'bowtie/3.4.5',
            version = '3.4.5',
        )
        Tool.objects.create(
            tool_name='Bowtie',
            tool_path = 'bowtie/2.4.5',
            version = '2.4.5',
        )
        Tool.objects.create(
            tool_name='Star',
            tool_path = 'star/3.4.5',
            version = '3.4.5',
        )
        res = Tool.objects.tools_list()
        assert res == {'Bowtie': ['2.4.5', '3.4.5'], 'Star': ['3.4.5']}

    def test_get_last_version(self):
        res = Tool.objects.get_last_version('abc')
        assert res is None

        Tool.objects.create(
            tool_name='Bowtie',
            tool_path = 'bowtie/2.0',
            version = '2.0',
        )
        Tool.objects.create(
            tool_name='Bowtie',
            tool_path = 'bowtie/1.0',
            version = '1.1.0',
        )
        Tool.objects.create(
            tool_name='Star',
            tool_path = 'star/3.4.5',
            version = '3.4.5',
        )
        res = Tool.objects.get_last_version('Bowtie')
        assert res.version == '2.0'

    def test_get_tool(self):
        res = Tool.objects.get_tool('abc', '1.0')
        assert res is None

        Tool.objects.create(
            tool_name='Bowtie',
            tool_path = 'bowtie/2.0',
            version = '2.0',
        )
        res = Tool.objects.get_tool('Bowtie', '2.0')
        assert res.version == '2.0'
        res = Tool.objects.get_tool('Bowtie', '6.0')
        assert res is None

    @override_settings(TOOLS_DIR='abc')
    def test_tools_path(self):
        assert Tool().working_dir == 'abc'

    @override_settings(TOOLS_DIR='/mnt/tools/')
    def test_full_too_path(self):
        t = Tool()
        t.tool_path = 'bowtie/3.4/bowtie'
        res = t.full_tool_path
        assert res == '/mnt/tools/bowtie/3.4/bowtie'
