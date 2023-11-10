import json
from django.db import models

class SampleManager(models.Manager):

    def study_exists(self, study_name:str)->bool:
        study = self.filter(study_name=study_name)
        if len(study) == 0:
            return False
        return True

    def sample_exists(self, study_name:str, sample_name:str)->bool:
        sample = self.filter(study_name=study_name, sample_name=sample_name)
        if len(sample) == 0:
            return False
        return True

    def get_sample_names_by_study(self, study_name:str)->list:
        '''
        return sample names given study_name
        '''
        samples = self.filter(study_name=study_name)
        return [s.sample_name for s in samples]

    def get_sample_names_by_user(self, user_name:str)->dict:
        '''
        given a user, return all studyes with sample names created by the user
        '''
        studyes = {}
        user = CustomUser.objects.get(user_name=user_name)
        samples = self.model.objects.filter(creator=user)
        for sample in samples:
            study_name = sample.study_name
            sample_name = sample.sample_name
            if study_name in studyes:
                studyes[study_name].append(sample_name)
            else:
                studyes[study_name] = [sample_name,]
        return studyes

    def get_study_names(self)->set:
        '''
        return all study names
        '''
        names = [i.study_name for i in self.all()]
        return set(names)
    
    def group_by_study(self)->set:
        '''
        group data by study name
        '''
        res = {}
        for i in self.all():
            if i.study_name in res:
                res[i.study_name].append(i)
            else:
                res[i.study_name] = [i,]
        return res
        
    def load_samples(self, user:str, samples:list):
        '''
        import samples into database
        study_name + sample_name should be unique
        '''
        res = {
            'created': [],
            'updated': [],
        }
        for sample in samples:
            existing_sample = self.filter(study_name = sample['study_name'], \
                sample_name=sample['sample_name'])
            print(existing_sample)
            if not existing_sample:
                new_sample = self.create(study_name = sample['study_name'], \
                    sample_name=sample['sample_name'], creator=user)
                if 'metadata' in sample:
                    new_sample.metadata = json.dumps(sample['metadata'])
                new_sample.save()
                res['created'].append(new_sample.id)
            else:
                if 'metadata' in sample:
                    existing_sample.update(metadata=json.dumps(sample['metadata']))
                    res['updated'] += [s.id for s in existing_sample]
        return res
    

    def update_sample_name(self, study_name:str, old_name:str, new_name:str):
        '''
        update sample name and keep unique
        '''
        sample = self.filter(study_name=study_name, sample_name=old_name)
        if len(sample) == 1:
            new_sample = self.model.objects.filter(
                study_name=study_name, sample_name=new_name)
            if len(new_sample) == 0:
                sample.update(sample_name=new_name)
                return self.get(study_name=study_name, sample_name=new_name)
        # elif len(sample) == 0:
        #     raise ValueError("The sample is not existing. " + \
        #         f"study_name={study_name}, sample_name={old_name}")
        return None

    def delete_sample(self, study_name:str, sample_name:str):
        return self.filter(study_name=study_name, sample_name=sample_name).delete()

    def delete_study_samples(self, study_name:str):
        return self.model.objects.filter(study_name=study_name).delete()

    def export_study(self, study_name:str)->dict:
        '''
        export table to nested dict
        '''
        study = {}
        samples = self.filter(study_name=study_name)
        for sample in samples:
            study[sample.sample_name]=vars(sample)
        return study


class Sample(models.Model):
    # one study include many samples
    study_name = models.CharField(max_length=50)
    # one sample on one name
    sample_name = models.CharField(max_length=100)
    creator = models.ForeignKey(
        'commons.CustomUser',
        on_delete=models.CASCADE
    )
    # namely phenotype
    # string type converted from json format
    metadata = models.CharField(max_length=3000, \
        blank=True, null=True)

    objects = SampleManager()
    unique_together = ('study_name', 'sample_name')

    class Meta:
        app_label = 'rna_seq'
        ordering = ('study_name', 'sample_name')

    def __str__(self):
        return f"{self.study_name}_{self.sample_name}"

    def to_dict(self):
        return {
            'study_name': self.study_name,
            'sample_name': self.sample_name,
            'creator': self.creator,
            'metadata': json.loads(self.metadata),
        }


