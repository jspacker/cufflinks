/*
 *  differential.cpp
 *  cufflinks
 *
 *  Created by Cole Trapnell on 3/15/10.
 *  Copyright 2010 Cole Trapnell. All rights reserved.
 *
 */

#include <algorithm>
#include <functional>
#include <numeric>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "abundances.h"
#include "differential.h"
#include "clustering.h"
#include "differential.h"
#include "sampling.h"

using namespace std;

double min_read_count = 10;
double min_outlier_p = 0.01;

#if ENABLE_THREADS
mutex _launcher_lock;
mutex locus_thread_pool_lock;
int locus_curr_threads = 0;
int locus_num_threads = 0;

void decr_pool_count()
{
	locus_thread_pool_lock.lock();
	locus_curr_threads--;
	locus_thread_pool_lock.unlock();	
}
#endif

void add_to_tracking_table(size_t sample_index,
                           Abundance& ab,
						   FPKMTrackingTable& track)

{
	pair<FPKMTrackingTable::iterator,bool> inserted;
	pair<string, FPKMTracking > p;
	p = make_pair(ab.description(), FPKMTracking());
	inserted = track.insert(p);
	
	FPKMTracking& fpkm_track = inserted.first->second;
	
	set<string> tss = ab.tss_id();
    set<string> gene_ids = ab.gene_id();
	set<string> genes = ab.gene_name();
	set<string> proteins = ab.protein_id();
	
	fpkm_track.tss_ids.insert(tss.begin(), tss.end());
    fpkm_track.gene_ids.insert(gene_ids.begin(), gene_ids.end());
	fpkm_track.gene_names.insert(genes.begin(), genes.end());
	fpkm_track.protein_ids.insert(proteins.begin(), proteins.end());
	
	if (inserted.second)
	{
		fpkm_track.locus_tag = ab.locus_tag();
		fpkm_track.description = ab.description();
		shared_ptr<Scaffold> transfrag = ab.transfrag();
		if (transfrag && transfrag->nearest_ref_id() != "")
		{
			fpkm_track.classcode = transfrag->nearest_ref_classcode();
			fpkm_track.ref_match = transfrag->nearest_ref_id();
		}
		else
		{
			fpkm_track.classcode = 0;
			fpkm_track.ref_match = "-";
		}
        if (transfrag)
        {
            fpkm_track.length = transfrag->length(); 
        }
        else
        {
            fpkm_track.length = 0;
        }
	}
	
	FPKMContext r1 = FPKMContext(ab.num_fragments(), 
                                 ab.num_fragment_var(),
                                 ab.num_fragment_uncertainty_var(),
                                 ab.mass_variance(),
                                 ab.num_fragments_by_replicate(),
								 ab.FPKM(), 
								 ab.FPKM_by_replicate(),
                                 ab.FPKM_variance(),
                                 ab.status(),
                                 ab.status_by_replicate());
    
    
	
    vector<FPKMContext>& fpkms = inserted.first->second.fpkm_series;
    if (sample_index < fpkms.size())
    {
        // if the fpkm series already has an entry matching this description
        // for this sample index, then we are dealing with a group of transcripts
        // that occupies multiple (genomically disjoint) bundles.  We need
        // to add this bundle's contribution to the FPKM, fragments, and variance 
        // to whatever's already there.  
        
        // NOTE: we can simply sum the FKPM_variances, because we are currently
        // assuming that transcripts in disjoint bundles share no alignments and 
        // thus have FPKM covariance == 0;  This assumption will no longer be
        // true if we decide to do multireads the right way.
        
        FPKMContext& existing = fpkms[sample_index];
        existing.FPKM += r1.FPKM;
        existing.count_mean += r1.count_mean;
        existing.FPKM_variance += r1.FPKM_variance;
        if (existing.status == NUMERIC_FAIL || r1.status == NUMERIC_FAIL)
        {
            existing.status = NUMERIC_FAIL;
        }
        else 
        {
            existing.status = NUMERIC_OK;
        }
        
    }
    else 
    {
        fpkms.push_back(r1);
    }
}

void add_to_allele_tracking_table(size_t sample_index,
                           Abundance& ab,
						   FPKMTrackingTable& track)

{
	pair<FPKMTrackingTable::iterator,bool> inserted;
	pair<string, FPKMTracking > p;
	p = make_pair(ab.description(), FPKMTracking());
	inserted = track.insert(p);
	
	FPKMTracking& fpkm_track = inserted.first->second;
	
	set<string> tss = ab.tss_id();
    set<string> gene_ids = ab.gene_id();
	set<string> genes = ab.gene_name();
	set<string> proteins = ab.protein_id();
	
	fpkm_track.tss_ids.insert(tss.begin(), tss.end());
    fpkm_track.gene_ids.insert(gene_ids.begin(), gene_ids.end());
	fpkm_track.gene_names.insert(genes.begin(), genes.end());
	fpkm_track.protein_ids.insert(proteins.begin(), proteins.end());
	
	if (inserted.second)
	{
		fpkm_track.locus_tag = ab.locus_tag();
		fpkm_track.description = ab.description();
		shared_ptr<Scaffold> transfrag = ab.transfrag();
		if (transfrag && transfrag->nearest_ref_id() != "")
		{
			fpkm_track.classcode = transfrag->nearest_ref_classcode();
			fpkm_track.ref_match = transfrag->nearest_ref_id();
		}
		else
		{
			fpkm_track.classcode = 0;
			fpkm_track.ref_match = "-";
		}
        if (transfrag)
        {
            fpkm_track.length = transfrag->length(); 
        }
        else
        {
            fpkm_track.length = 0;
        }
	}
	

	FPKMContext paternal_r1 = FPKMContext(ab.num_paternal_fragments(), 
                                 ab.num_paternal_fragment_var(),
                                 ab.num_paternal_fragment_uncertainty_var(),
                                 ab.paternal_mass_variance(),
                                 ab.num_paternal_fragments_by_replicate(),
								 ab.paternal_FPKM(), 
								 ab.paternal_FPKM_by_replicate(),
                                 ab.paternal_FPKM_variance(),
                                 ab.paternal_status(),
                                 ab.paternal_status_by_replicate());
	FPKMContext maternal_r1 = FPKMContext(ab.num_maternal_fragments(), 
                                 ab.num_maternal_fragment_var(),
                                 ab.num_maternal_fragment_uncertainty_var(),
                                 ab.maternal_mass_variance(),
                                 ab.num_maternal_fragments_by_replicate(),
								 ab.maternal_FPKM(), 
								 ab.maternal_FPKM_by_replicate(),
                                 ab.maternal_FPKM_variance(),
                                 ab.maternal_status(),
                                 ab.maternal_status_by_replicate());
    
	vector<FPKMContext>& paternal_fpkms = inserted.first->second.paternal_fpkm_series;
    if (sample_index < paternal_fpkms.size())
    {
        // if the fpkm series already has an entry matching this description
        // for this sample index, then we are dealing with a group of transcripts
        // that occupies multiple (genomically disjoint) bundles.  We need
        // to add this bundle's contribution to the FPKM, fragments, and variance 
        // to whatever's already there.  
        
        // NOTE: we can simply sum the FKPM_variances, because we are currently
        // assuming that transcripts in disjoint bundles share no alignments and 
        // thus have FPKM covariance == 0;  This assumption will no longer be
        // true if we decide to do multireads the right way.
        
        FPKMContext& existing = paternal_fpkms[sample_index];
		
        existing.FPKM += paternal_r1.FPKM;
        existing.count_mean += paternal_r1.count_mean;
        existing.FPKM_variance += paternal_r1.FPKM_variance;
        if (existing.status == NUMERIC_FAIL || paternal_r1.status == NUMERIC_FAIL)
        {
            existing.status = NUMERIC_FAIL;
        }
        else 
        {
            existing.status = NUMERIC_OK;
        }
        
    }
    else 
    {
        paternal_fpkms.push_back(paternal_r1);
    }
	vector<FPKMContext>& maternal_fpkms = inserted.first->second.maternal_fpkm_series;
    if (sample_index < maternal_fpkms.size())
    {
        // if the fpkm series already has an entry matching this description
        // for this sample index, then we are dealing with a group of transcripts
        // that occupies multiple (genomically disjoint) bundles.  We need
        // to add this bundle's contribution to the FPKM, fragments, and variance 
        // to whatever's already there.  
        
        // NOTE: we can simply sum the FKPM_variances, because we are currently
        // assuming that transcripts in disjoint bundles share no alignments and 
        // thus have FPKM covariance == 0;  This assumption will no longer be
        // true if we decide to do multireads the right way.
        
        FPKMContext& existing = maternal_fpkms[sample_index];
		
        existing.FPKM += maternal_r1.FPKM;
        existing.count_mean += maternal_r1.count_mean;
        existing.FPKM_variance += maternal_r1.FPKM_variance;
        if (existing.status == NUMERIC_FAIL || maternal_r1.status == NUMERIC_FAIL)
        {
            existing.status = NUMERIC_FAIL;
        }
        else 
        {
            existing.status = NUMERIC_OK;
        }
        
    }
    else 
    {
        maternal_fpkms.push_back(maternal_r1);
    }
}


TestLauncher::launcher_sample_table::iterator TestLauncher::find_locus(const string& locus_id)
{
    launcher_sample_table::iterator itr = _samples.begin();
    for(; itr != _samples.end(); ++itr)
    {
        if (itr->first == locus_id)
            return itr;
    }
    return _samples.end();
}

void TestLauncher::register_locus(const string& locus_id)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif	
    
    launcher_sample_table::iterator itr = find_locus(locus_id);
    if (itr == _samples.end())
    {
        pair<launcher_sample_table::iterator, bool> p;
        vector<shared_ptr<SampleAbundances> >abs(_orig_workers);
        _samples.push_back(make_pair(locus_id, abs));
    }
}

void TestLauncher::abundance_avail(const string& locus_id, 
                                   shared_ptr<SampleAbundances> ab, 
                                   size_t factory_id)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif	
    launcher_sample_table::iterator itr = find_locus(locus_id);
    if (itr == _samples.end())
    {
        assert(false);
    }
    itr->second[factory_id] = ab;
    //itr->second(factory_id] = ab;
}

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
bool TestLauncher::all_samples_reported_in(vector<shared_ptr<SampleAbundances> >& abundances)
{    
    foreach (shared_ptr<SampleAbundances> ab, abundances)
    {
        if (!ab)
        {
            return false;
        }
    }
    return true;
}

#if ENABLE_THREADS
mutex test_storage_lock; // don't modify the above struct without locking here
#endif

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
void TestLauncher::perform_testing(vector<shared_ptr<SampleAbundances> >& abundances)
{
    assert (abundances.size() == _orig_workers);
    
    // Just verify that all the loci from each factory match up.
    for (size_t i = 1; i < abundances.size(); ++i)
    {
        const SampleAbundances& curr = *(abundances[i]);
        const SampleAbundances& prev = *(abundances[i-1]);
        
        assert (curr.locus_tag == prev.locus_tag);
        
        const AbundanceGroup& s1 = curr.transcripts;
        const AbundanceGroup& s2 =  prev.transcripts;
        
        assert (s1.abundances().size() == s2.abundances().size());
        
        for (size_t j = 0; j < s1.abundances().size(); ++j)
        {
            assert (s1.abundances()[j]->description() == s2.abundances()[j]->description());
        }
    }
    
    test_differential(abundances.front()->locus_tag, abundances, *_tests, *_tracking, _samples_are_time_series);
}

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
void TestLauncher::record_tracking_data(vector<shared_ptr<SampleAbundances> >& abundances)
{
    assert (abundances.size() == _orig_workers);
    
    // Just verify that all the loci from each factory match up.
    for (size_t i = 1; i < abundances.size(); ++i)
    {
        const SampleAbundances& curr = *(abundances[i]);
        const SampleAbundances& prev = *(abundances[i-1]);
        
        assert (curr.locus_tag == prev.locus_tag);
        
        const AbundanceGroup& s1 = curr.transcripts;
        const AbundanceGroup& s2 =  prev.transcripts;
        
        assert (s1.abundances().size() == s2.abundances().size());
        
        for (size_t j = 0; j < s1.abundances().size(); ++j)
        {
            assert (s1.abundances()[j]->description() == s2.abundances()[j]->description());
        }
    }
    
#if ENABLE_THREADS
	test_storage_lock.lock();
#endif
    
    // Add all the transcripts, CDS groups, TSS groups, and genes to their
    // respective FPKM tracking table.  Whether this is a time series or an
    // all pairs comparison, we should be calculating and reporting FPKMs for 
    // all objects in all samples
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		const AbundanceGroup& ab_group = abundances[i]->transcripts;
        //fprintf(stderr, "[%d] count = %lg\n",i,  ab_group.num_fragments());
		foreach (shared_ptr<Abundance> ab, ab_group.abundances())
		{
			add_to_tracking_table(i, *ab, _tracking->isoform_fpkm_tracking);
            //assert (_tracking->isoform_fpkm_tracking.num_fragments_by_replicate().empty() == false);
		}
		
		foreach (AbundanceGroup& ab, abundances[i]->cds)
		{
			add_to_tracking_table(i, ab, _tracking->cds_fpkm_tracking);
		}
		
		foreach (AbundanceGroup& ab, abundances[i]->primary_transcripts)
		{
			add_to_tracking_table(i, ab, _tracking->tss_group_fpkm_tracking);
		}
		
		foreach (AbundanceGroup& ab, abundances[i]->genes)
		{
			add_to_tracking_table(i, ab, _tracking->gene_fpkm_tracking);
		}
	}
	
	
#if ENABLE_THREADS
    test_storage_lock.unlock();
#endif
    
}

void TestLauncher::test_finished_loci()
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif  

    launcher_sample_table::iterator itr = _samples.begin(); 
    while(itr != _samples.end())
    {
        if (all_samples_reported_in(itr->second))
        {
            // In some abundance runs, we don't actually want to perform testing 
            // (eg initial quantification before bias correction).
            // _tests and _tracking will be NULL in these cases.
            if (_tests != NULL && _tracking != NULL)
            {
                if (_p_bar)
                {
                    verbose_msg("Testing for differential expression and regulation in locus [%s]\n", itr->second.front()->locus_tag.c_str());
					_p_bar->update(itr->second.front()->locus_tag.c_str(), 1);
				}
                record_tracking_data(itr->second);
				perform_testing(itr->second);
            }
            else
            {
                if (_p_bar)
                {
                    //verbose_msg("Testing for differential expression and regulation in locus [%s]\n", abundances.front()->locus_tag.c_str());
                    _p_bar->update(itr->second.front()->locus_tag.c_str(), 1);
                }
                record_tracking_data(itr->second);
            }
            itr = _samples.erase(itr);
        }
        else
        {
            
            ++itr;
        }
    }
}


AlleleTestLauncher::launcher_sample_table::iterator AlleleTestLauncher::find_locus(const string& locus_id)
{
    launcher_sample_table::iterator itr = _samples.begin();
    for(; itr != _samples.end(); ++itr)
    {
        if (itr->first == locus_id)
            return itr;
    }
    return _samples.end();
}

void AlleleTestLauncher::register_locus(const string& locus_id)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif	
    
    launcher_sample_table::iterator itr = find_locus(locus_id);
    if (itr == _samples.end())
    {
        pair<launcher_sample_table::iterator, bool> p;
        vector<shared_ptr<SampleAlleleAbundances> >abs(_orig_workers);
        _samples.push_back(make_pair(locus_id, abs));
    }
}

void AlleleTestLauncher::abundance_avail(const string& locus_id, 
                                   shared_ptr<SampleAlleleAbundances> ab, 
                                   size_t factory_id)
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif	
    launcher_sample_table::iterator itr = find_locus(locus_id);
    if (itr == _samples.end())
    {
        assert(false);
    }
    itr->second[factory_id] = ab;
    //itr->second(factory_id] = ab;
}

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
bool AlleleTestLauncher::all_samples_reported_in(vector<shared_ptr<SampleAlleleAbundances> >& abundances)
{    
    foreach (shared_ptr<SampleAlleleAbundances> ab, abundances)
    {
        if (!ab)
        {
            return false;
        }
    }
    return true;
}

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
void AlleleTestLauncher::perform_testing(vector<shared_ptr<SampleAlleleAbundances> >& abundances)
{
    assert (abundances.size() == _orig_workers);
	// Just verify that all the loci from each factory match up.
	if(abundances.size() > 1) //this means that we are testing each combination of alleles in each sample/condition (not the classic imprinting case in which we sinply have 1 sample/condition)
	{
		for (size_t i = 1; i < abundances.size(); ++i)
		{
			const SampleAlleleAbundances& curr = *(abundances[i]);
			const SampleAlleleAbundances& prev = *(abundances[i-1]);
			
			assert (curr.locus_tag == prev.locus_tag);
			
			const AlleleAbundanceGroup& s1 = curr.transcripts;
			const AlleleAbundanceGroup& s2 =  prev.transcripts;
			
			assert (s1.abundances().size() == s2.abundances().size());
			
			for (size_t j = 0; j < s1.abundances().size(); ++j)
			{
				assert (s1.abundances()[j]->description() == s2.abundances()[j]->description());
			}
		}
	}
	test_differential(abundances.front()->locus_tag, abundances, *_tests, *_tracking, _samples_are_time_series);
}

// Note: this routine should be called under lock - it doesn't
// acquire the lock itself. 
void AlleleTestLauncher::record_tracking_data(vector<shared_ptr<SampleAlleleAbundances> >& abundances)
{
    assert (abundances.size() == _orig_workers);
	// Just verify that all the loci from each factory match up.
    for (size_t i = 1; i < abundances.size(); ++i)
    {
        const SampleAlleleAbundances& curr = *(abundances[i]);
        const SampleAlleleAbundances& prev = *(abundances[i-1]);
        
        assert (curr.locus_tag == prev.locus_tag);
        
        const AlleleAbundanceGroup& s1 = curr.transcripts;
        const AlleleAbundanceGroup& s2 = prev.transcripts;
        
        assert (s1.abundances().size() == s2.abundances().size());
        
        for (size_t j = 0; j < s1.abundances().size(); ++j)
        {
            assert (s1.abundances()[j]->description() == s2.abundances()[j]->description());
        }
    }
#if ENABLE_THREADS
	test_storage_lock.lock();
#endif
    
    // Add all the transcripts, CDS groups, TSS groups, and genes to their
    // respective FPKM tracking table.  Whether this is a time series or an
    // all pairs comparison, we should be calculating and reporting FPKMs for 
    // all objects in all samples
	for (size_t i = 0; i < abundances.size(); ++i)
	{
		const AlleleAbundanceGroup& ab_group = abundances[i]->transcripts;
        //fprintf(stderr, "[%d] count = %lg\n",i,  ab_group.num_fragments());
		//nimrod??
		foreach (shared_ptr<Abundance> ab, ab_group.abundances())
		{
			add_to_allele_tracking_table(i, *ab, _tracking->isoform_fpkm_tracking);
            //assert (_tracking->isoform_fpkm_tracking.num_fragments_by_replicate().empty() == false);
		}
		foreach (AlleleAbundanceGroup& ab, abundances[i]->cds)
		{
			add_to_allele_tracking_table(i, ab, _tracking->cds_fpkm_tracking);
		}
		foreach (AlleleAbundanceGroup& ab, abundances[i]->primary_transcripts)
		{
			add_to_allele_tracking_table(i, ab, _tracking->tss_group_fpkm_tracking);
		}
		foreach (AlleleAbundanceGroup& ab, abundances[i]->genes)
		{
			add_to_allele_tracking_table(i, ab, _tracking->gene_fpkm_tracking);
		}
	}
#if ENABLE_THREADS
    test_storage_lock.unlock();
#endif
}

void AlleleTestLauncher::test_finished_loci()
{
#if ENABLE_THREADS
	boost::mutex::scoped_lock lock(_launcher_lock);
#endif  

    launcher_sample_table::iterator itr = _samples.begin(); 
    while(itr != _samples.end())
    {
		if (all_samples_reported_in(itr->second))
        {
            // In some abundance runs, we don't actually want to perform testing 
            // (eg initial quantification before bias correction).
            // _tests and _tracking will be NULL in these cases.
            if (_tests != NULL && _tracking != NULL)
            {
                if (_p_bar)
                {
                    verbose_msg("Testing for allele differential expression and regulation in locus [%s]\n", itr->second.front()->locus_tag.c_str());
					_p_bar->update(itr->second.front()->locus_tag.c_str(), 1);
				}
				record_tracking_data(itr->second);
				perform_testing(itr->second);
			}
            else
            {
                if (_p_bar)
                {
                    //verbose_msg("Testing for allele differential expression and regulation in locus [%s]\n", abundances.front()->locus_tag.c_str());
                    _p_bar->update(itr->second.front()->locus_tag.c_str(), 1);
                }
                record_tracking_data(itr->second);
            }
            itr = _samples.erase(itr);
        }
        else
        {
            
            ++itr;
        }
    }
}

// This performs a between-group test on an isoform or TSS grouping, on two 
// different samples.
bool test_diffexp(const FPKMContext& curr,
				  const FPKMContext& prev,
				  SampleDifference& test)
{
	bool performed_test = false;
	if (curr.FPKM > 0.0 && prev.FPKM > 0.0)
	{
		//assert (curr.FPKM_variance > 0.0 && prev.FPKM_variance > 0.0);
//		double log_curr = log(curr.counts);
//		double log_prev = log(prev.counts);
        
        double stat = 0.0;
        double p_value = 1.0;
				
        if (curr.FPKM_variance > 0.0 || prev.FPKM_variance > 0.0)
        {
            double curr_log_fpkm_var = (curr.FPKM_variance) / (curr.FPKM * curr.FPKM);
            double prev_log_fpkm_var = (prev.FPKM_variance) / (prev.FPKM * prev.FPKM);
            
            double numerator = log(prev.FPKM / curr.FPKM);
            
            double denominator = sqrt(prev_log_fpkm_var + curr_log_fpkm_var);
            stat = numerator / denominator;
        
		
            normal norm;
            double t1, t2;
            if (stat > 0.0)
            {
                t1 = stat;
                t2 = -stat;
            }
            else
            {
                t1 = -stat;
                t2 = stat;
            }
			if (isnan(t1) || isinf(t1) || isnan(t2) || isnan(t2))
            {
                
                //fprintf(stderr, "Warning: test statistic is NaN! %s (samples %lu and %lu)\n", test.locus_desc.c_str(), test.sample_1, test.sample_2);
				p_value = 1.0;
            }
            else
            {
                double tail_1 = cdf(norm, t1);
                double tail_2 = cdf(norm, t2);
                p_value = 1.0 - (tail_1 - tail_2);                
            }
        }
		
		double differential = log2(curr.FPKM) - log2(prev.FPKM);
		
		//test = SampleDifference(sample1, sample2, prev.FPKM, curr.FPKM, stat, p_value, transcript_group_id);
		test.p_value = p_value;
		test.differential = differential;
		test.test_stat = stat;
		test.value_1 = prev.FPKM;
		test.value_2 = curr.FPKM;
		
		performed_test = true;
	}
	else
	{
		if (curr.FPKM > 0.0)
		{
            if (curr.status != NUMERIC_LOW_DATA && curr.FPKM_variance > 0.0)
            {
                normal norm(curr.FPKM, sqrt(curr.FPKM_variance));
                test.p_value = cdf(norm, 0);
                performed_test = true;
                test.differential = numeric_limits<double>::max();;
                test.test_stat = numeric_limits<double>::max();
                test.value_1 = 0;
                test.value_2 = curr.FPKM;
            }
            else
            {
                test.differential = -numeric_limits<double>::max();
                test.test_stat = -numeric_limits<double>::max();
                test.value_1 = prev.FPKM;
                test.value_2 = 0;
                test.p_value = 1;
                performed_test = false;
            }
		}
		else if (prev.FPKM > 0.0)
		{
            if (prev.status != NUMERIC_LOW_DATA && prev.FPKM_variance > 0.0)
            {
                normal norm(prev.FPKM, sqrt(prev.FPKM_variance));
                test.p_value = cdf(norm, 0);
                performed_test = true;
                
                test.differential = -numeric_limits<double>::max();
                test.test_stat = -numeric_limits<double>::max();
                test.value_1 = prev.FPKM;
                test.value_2 = 0;
            }
            else
            {
                test.differential = -numeric_limits<double>::max();
                test.test_stat = -numeric_limits<double>::max();
                test.value_1 = prev.FPKM;
                test.value_2 = 0;
                test.p_value = 1;
                performed_test = false;
            }
		}
        else
        {
            assert (prev.FPKM == 0.0 && curr.FPKM == 0.0);
            performed_test = false;
        }
        
	}	
	
	test.test_status = performed_test ? OK : NOTEST;
	return performed_test;
}

SampleDiffMetaDataTable meta_data_table;
#if ENABLE_THREADS
boost::mutex meta_data_lock;
#endif

shared_ptr<SampleDifferenceMetaData> get_metadata(const string description)
{
#if ENABLE_THREADS
    boost::mutex::scoped_lock lock(meta_data_lock);
#endif
    pair<SampleDiffMetaDataTable::iterator, bool> p;
    p = meta_data_table.insert(make_pair(description, new SampleDifferenceMetaData()));
    return p.first->second;
}

// This performs between-group tests on isoforms or TSS groupings in a single
// locus, on two different samples.
SampleDifference get_de_tests(const string& description,
                 const FPKMContext& prev_abundance,
				 const FPKMContext& curr_abundance,
				 //SampleDiffs& de_tests,
							  bool enough_reads,
							  bool doTest=true)
{
	int total_iso_de_tests = 0;
			
	SampleDifference test;
    
    const FPKMContext& r1 = curr_abundance;
    const FPKMContext& r2 = prev_abundance;
    
	if(!doTest)
	{
		test.test_stat = 0;
		test.p_value = 1.0;
		test.differential = 0.0;
		test.test_status = NOALLELETEST;
		return test;
	}
	if (curr_abundance.status == NUMERIC_FAIL || 
        prev_abundance.status == NUMERIC_FAIL ||
        prev_abundance.status == NUMERIC_HI_DATA ||
        curr_abundance.status == NUMERIC_HI_DATA)
    {
		test_diffexp(r1, r2, test);
        test.test_stat = 0;
		test.p_value = 1.0;
		test.differential = 0.0;
        if (curr_abundance.status == NUMERIC_FAIL || 
            prev_abundance.status == NUMERIC_FAIL)
        {
            test.test_status = FAIL;
        }
        else if (prev_abundance.status == NUMERIC_HI_DATA ||
                 curr_abundance.status == NUMERIC_HI_DATA)
        {
            test.test_status = HIDATA;
        }
    }
    else if (curr_abundance.status == NUMERIC_LOW_DATA || 
             prev_abundance.status == NUMERIC_LOW_DATA)
    {
		// perform the test, but mark it as not significant and don't add it to the 
        // pile. This way we don't penalize for multiple testing, but users can still
        // see the fold change.
		test_diffexp(r1, r2, test);
        test.test_stat = 0;
        test.p_value = 1.0;
        //test.differential = 0.0;
		
		test.test_status = LOWDATA;
    }
    else // at least one is OK, the other might be LOW_DATA
	{
		test.test_status = FAIL;
		if (test_diffexp(r1, r2, test))
		{
			total_iso_de_tests++;
		}
		else
		{
			test.test_stat = 0;
			test.p_value = 1.0;
			test.differential = 0.0;
		}
		if (enough_reads) 
			test.test_status = OK;
		else
			test.test_status = NOTEST;
		
	}
	
	//inserted.first->second = test;
	
	//return make_pair(total_iso_de_tests, inserted.first);
	return test;
}



//bool generate_js_samples(const AbundanceGroup& prev_abundance,
//                         const AbundanceGroup& curr_abundance,
//                         size_t num_js_samples,
//                         vector<double>& js_samples)
//{
//    ublas::vector<double> prev_kappa_mean(curr_abundance.abundances().size());
//    for (size_t i = 0; i < prev_abundance.abundances().size(); ++i)
//    {
//        prev_kappa_mean(i) = prev_abundance.abundances()[i]->kappa();
//    }
//    
//    ublas::vector<double> curr_kappa_mean(curr_abundance.abundances().size());
//    for (size_t i = 0; i < curr_abundance.abundances().size(); ++i)
//    {
//        curr_kappa_mean(i) = curr_abundance.abundances()[i]->kappa();
//    }
//
//    ublas::matrix<double> prev_kappa_cov = prev_abundance.kappa_cov();
//    double prev_ret = cholesky_factorize(prev_kappa_cov);
//    if (prev_ret != 0)
//        return false;
//    
//    ublas::matrix<double> curr_kappa_cov = curr_abundance.kappa_cov();
//    double curr_ret = cholesky_factorize(curr_kappa_cov);
//    if (curr_ret != 0)
//        return false;
//    
//    vector<ublas::vector<double> > samples;
//    
//    multinormal_generator<double> generator(prev_kappa_mean, prev_kappa_cov);
//    //vector<ublas::vector<double> > prev_samples;
//    generate_importance_samples(generator, samples, num_js_samples / 2, true);
//    
//    // It's a little silly that we have to do this, but since we always initialize
//    // the random number generators to random_seed, instead of time(NULL), simply
//    // creating a new generator (rather than re-using it) 
//    generator.set_parameters(curr_kappa_mean, curr_kappa_cov);
//    
//    //multinormal_generator<double> curr_generator(curr_kappa_mean, curr_kappa_cov);
//    //vector<ublas::vector<double> > curr_samples;
//    generate_importance_samples(generator, samples, num_js_samples / 2, true);
//
//    // We want to revise the covariance matrix from the samples, since we'll 
//    // need it later for the CIs.
//    ublas::matrix<double> null_kappa_cov;
//    null_kappa_cov = ublas::zero_matrix<double>(curr_kappa_cov.size1(),
//                                               curr_kappa_cov.size2());
//    
//
//    ublas::vector<double> null_kappa_mean;
//    null_kappa_mean = ublas::zero_vector<double>(curr_kappa_cov.size1());
//    
//    foreach(ublas::vector<double>& sample, samples)
//    {
//        null_kappa_mean += sample;
//    }
//    null_kappa_mean /= samples.size();
//   
//    for (size_t i = 0; i < null_kappa_cov.size1(); ++i)
//    {
//        for (size_t j = 0; j < null_kappa_cov.size2(); ++j)
//        {
//            for (size_t k = 0 ; k < samples.size(); ++k)
//            {
//                double c = (samples[k](i) - null_kappa_mean(i)) * (samples[k](j) - null_kappa_mean(j));
//                null_kappa_cov(i,j) += c;
//            }
//        }
//    }
//    
//    null_kappa_cov /= samples.size();
//    
//    static const double epsilon = 1e-6;
//    
//    null_kappa_cov += (epsilon * ublas::identity_matrix<double>(null_kappa_cov.size1()));
//    
//    double null_ret = cholesky_factorize(null_kappa_cov);
//    if (null_ret != 0)
//        return false;
//    
//    generator.set_parameters(null_kappa_mean, null_kappa_cov);
//    samples.clear();
//    generate_importance_samples(generator, samples, num_js_samples, true);
//    
//    cerr << "prev: " << endl;
//    cerr << prev_kappa_mean << endl;
//    for (unsigned i = 0; i < prev_kappa_cov.size1 (); ++ i) 
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (prev_kappa_cov, i);
//        std::cerr << i << " : " << mr << std::endl;
//    }
//    cerr << "curr: " << endl;
//    cerr << curr_kappa_mean << endl;
//    for (unsigned i = 0; i < curr_kappa_cov.size1 (); ++ i) 
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (curr_kappa_cov, i);
//        std::cerr << i << " : " << mr << std::endl;
//    }
//    cerr << "null: " << endl;
//    cerr << null_kappa_mean << endl;
//    for (unsigned i = 0; i < null_kappa_cov.size1 (); ++ i) 
//    {
//        ublas::matrix_row<ublas::matrix<double> > mr (null_kappa_cov, i);
//        std::cerr << i << " : " << mr << std::endl;
//    }
//    cerr << "======" << endl;
//    
//    js_samples.clear();
//    
//    
//    //size_t num_samples = std::min(prev_samples.size(), curr_samples.size());
//    size_t num_samples = num_js_samples;
//    vector<ublas::vector<double> > sample_kappas(2);
//    
//    boost::uniform_int<> uniform_dist(0,samples.size()-1);
//    boost::mt19937 rng; 
//    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > uniform_gen(rng, uniform_dist); 
//    
//    for (size_t i = 0; i < num_samples; ++i)
//    {
//		sample_kappas[0] = samples[uniform_gen()];
//        sample_kappas[1] = samples[uniform_gen()];
//            
//		double js = jensen_shannon_distance(sample_kappas);  
//        cerr << sample_kappas[0] << " vs. " <<  sample_kappas[1] << " = " << js << endl;
//        js_samples.push_back(js);
//    }
//    
//    sort(js_samples.begin(), js_samples.end());
//    
////    for (size_t i = 0; i < 100; ++i)
////    {
////        fprintf(stderr, "%lg\n", js_samples[i]);
////    }
//    return true;
//}

// Calculates the probability that drawing two samples from the provided
// relative abundance distribution would have produced a value at least as
// extreme as the given js value.
bool one_sided_js_test(const AbundanceGroup& null_abundances,
                       double js,
                       double& p_val)
{
    const vector<double>& js_samples = null_abundances.null_js_samples();
    
    vector<double>::const_iterator lb = lower_bound(js_samples.begin(), js_samples.end(), js);
    if (lb != js_samples.end())
    {
        size_t num_less_extreme_samples = lb - js_samples.begin();
        p_val =  1.0  - ((double)num_less_extreme_samples/js_samples.size());
    }
    else if (js_samples.size())
    {
        p_val = 1.0/js_samples.size();
    }
    else
    {
        p_val = 1.0;
        return false;
    }
    return true;
}

// Calculates the probability that drawing two samples from the provided
// relative abundance distribution would have produced a value at least as
// extreme as the given js value.
bool one_sided_js_test(const vector<double>& js_samples,
                       double js,
                       double& p_val)
{
	vector<double>::const_iterator lb = lower_bound(js_samples.begin(), js_samples.end(), js);
    if (lb != js_samples.end())
    {
        size_t num_less_extreme_samples = lb - js_samples.begin();
        p_val =  1.0  - ((double)num_less_extreme_samples/js_samples.size());
    }
    else if (js_samples.size())
    {
        p_val = 1.0/js_samples.size();
    }
    else
    {
        p_val = 1.0;
        return false;
    }
    return true;
}

bool test_js(const AbundanceGroup& prev_abundance,
             const AbundanceGroup& curr_abundance,
             double& js,
             double& p_val)
{
    vector<Eigen::VectorXd> sample_kappas;
    Eigen::VectorXd curr_kappas(Eigen::VectorXd::Zero(curr_abundance.abundances().size()));
    for (size_t i = 0; i < curr_abundance.abundances().size(); ++i)
    {
        curr_kappas(i) = curr_abundance.abundances()[i]->kappa();
    }
    
    Eigen::VectorXd prev_kappas(Eigen::VectorXd::Zero(prev_abundance.abundances().size()));
    for (size_t i = 0; i < prev_abundance.abundances().size(); ++i)
    {
        prev_kappas(i) = prev_abundance.abundances()[i]->kappa();
    }
    
    sample_kappas.push_back(prev_kappas);
    sample_kappas.push_back(curr_kappas);
    
    js = jensen_shannon_distance(sample_kappas);
    
    if (isinf(js) || isnan(js))
        return false;
    
    double prev_p_val = 1.0;
    bool prev_succ = one_sided_js_test(prev_abundance, js, prev_p_val);
    
    double curr_p_val = 1.0;
    bool curr_succ = one_sided_js_test(curr_abundance, js, curr_p_val);

    if (!curr_succ || !prev_succ)
        return false;
        
    p_val = (prev_p_val + curr_p_val)/2;
//    double mean_js = accumulate(js_samples.begin(), js_samples.end(), 0.0);
//    if (js_samples.size() == 0)
//        return false;
//    mean_js /= js_samples.size();
//    
//    double var_js = 0.0;
//    for (size_t i = 0; i < js_samples.size(); ++i)
//    {
//        double s =  js_samples[i] - mean_js;
//        s *= s;
//        var_js += s;
//    }
//    var_js /= js_samples.size();
//    
//    if (var_js > 0.0)
//    {
//        // We're dealing with a standard normal that's been truncated below zero
//        // so pdf(js) is twice the standard normal, and cdf is 0.5 * (cdf of normal - 1)
//        
//        normal test_dist(0,1.0);
//        //double denom = sqrt(js_var);
//        double se_js = sqrt(var_js);
//        double p = mean_js/se_js;
//        p_val = 1.0 - ((cdf(test_dist, p) - 0.5) / 0.5);
//    }
//    else
//    {
//        return false;
//    }
    
    return true;
}

bool test_js(const AlleleAbundanceGroup& prev_abundance,
             const AlleleAbundanceGroup& curr_abundance,
             double& js,
             double& p_val,
			 bool prev_parent,
			 bool curr_parent)
{
	//0-prev_parent,curr_parent = paternal,1-prev_parent,curr_parent = maternal
    vector<Eigen::VectorXd> sample_kappas;
    Eigen::VectorXd curr_kappas(Eigen::VectorXd::Zero(curr_abundance.abundances().size()));
    for (size_t i = 0; i < curr_abundance.abundances().size(); ++i)
    {
		if(!curr_parent)
		{
			curr_kappas(i) = curr_abundance.abundances()[i]->paternal_kappa();
		}
		else
		{
			curr_kappas(i) = curr_abundance.abundances()[i]->maternal_kappa();
		}
    }
    
    Eigen::VectorXd prev_kappas(Eigen::VectorXd::Zero(prev_abundance.abundances().size()));
    for (size_t i = 0; i < prev_abundance.abundances().size(); ++i)
    {
		if(!prev_parent)
		{
			prev_kappas(i) = prev_abundance.abundances()[i]->paternal_kappa();
		}
		else
		{
			prev_kappas(i) = prev_abundance.abundances()[i]->maternal_kappa();
		}
    }
    
    sample_kappas.push_back(prev_kappas);
    sample_kappas.push_back(curr_kappas);
    
    js = jensen_shannon_distance(sample_kappas);
    
    if (isinf(js) || isnan(js))
        return false;
    
    double prev_p_val = 1.0;
	bool prev_succ,curr_succ;
	if(!prev_parent)
		prev_succ = one_sided_js_test(prev_abundance.paternal_null_js_samples(), js, prev_p_val);
	else
		prev_succ = one_sided_js_test(prev_abundance.maternal_null_js_samples(), js, prev_p_val);
    double curr_p_val = 1.0;
	if(!curr_parent)
		curr_succ = one_sided_js_test(curr_abundance.paternal_null_js_samples(), js, curr_p_val);
	else
		curr_succ = one_sided_js_test(curr_abundance.maternal_null_js_samples(), js, curr_p_val);

    if (!curr_succ || !prev_succ)
        return false;
        
    p_val = (prev_p_val + curr_p_val)/2;
//    double mean_js = accumulate(js_samples.begin(), js_samples.end(), 0.0);
//    if (js_samples.size() == 0)
//        return false;
//    mean_js /= js_samples.size();
//    
//    double var_js = 0.0;
//    for (size_t i = 0; i < js_samples.size(); ++i)
//    {
//        double s =  js_samples[i] - mean_js;
//        s *= s;
//        var_js += s;
//    }
//    var_js /= js_samples.size();
//    
//    if (var_js > 0.0)
//    {
//        // We're dealing with a standard normal that's been truncated below zero
//        // so pdf(js) is twice the standard normal, and cdf is 0.5 * (cdf of normal - 1)
//        
//        normal test_dist(0,1.0);
//        //double denom = sqrt(js_var);
//        double se_js = sqrt(var_js);
//        double p = mean_js/se_js;
//        p_val = 1.0 - ((cdf(test_dist, p) - 0.5) / 0.5);
//    }
//    else
//    {
//        return false;
//    }
    
    return true;
}

// This performs within-group tests on a set of isoforms or a set of TSS groups.
// This is a way of looking for meaningful differential splicing or differential
// promoter use.
SampleDifference get_ds_tests(const AbundanceGroup& prev_abundance,
                              const AbundanceGroup& curr_abundance,
//                              SampleDiffs& diff_tests,
                              bool enough_reads)
{	
	const string& name = curr_abundance.description();
	
	SampleDifference test;
    
	test.test_status = NOTEST;
	
	AbundanceStatus prev_status = curr_abundance.status();
	AbundanceStatus curr_status = prev_abundance.status();
    
    if (prev_abundance.abundances().size() == 1 ||
        (prev_status == NUMERIC_OK && prev_abundance.num_fragments() == 0) ||
        (curr_status == NUMERIC_OK && curr_abundance.num_fragments() == 0))
    {
        test.p_value = 1;
        test.value_1 = 0;
        test.value_2 = 0;
        test.differential = 0;
        test.test_status = NOTEST;
    }
	else if (prev_abundance.abundances().size() > 1 &&
        /*prev_abundance.has_member_with_status(NUMERIC_LOW_DATA) == false &&
        filtered_curr.has_member_with_status(NUMERIC_LOW_DATA) == false &&*/ 
        prev_status == NUMERIC_OK && prev_abundance.num_fragments() > 0 &&
		curr_status == NUMERIC_OK && curr_abundance.num_fragments() > 0)
	{
		vector<ublas::vector<double> > sample_kappas;
		ublas::vector<double> curr_kappas(curr_abundance.abundances().size());
        
        for (size_t i = 0; i < curr_abundance.abundances().size(); ++i)
		{
			curr_kappas(i) = curr_abundance.abundances()[i]->kappa();
		}
		
		ublas::vector<double> prev_kappas(prev_abundance.abundances().size());
        for (size_t i = 0; i < prev_abundance.abundances().size(); ++i)
		{
			prev_kappas(i) = prev_abundance.abundances()[i]->kappa();
		}
		
		sample_kappas.push_back(prev_kappas);
		sample_kappas.push_back(curr_kappas);
		
        double js = 0.0;
        double p_val = 1.0;
        
        bool success;
        success = test_js(prev_abundance, curr_abundance, js, p_val);
        
		if (js == 0.0 || success == false)
		{
			test.test_stat = 0;
			test.p_value = 1.0;
			test.value_1 = 0;
			test.value_2 = 0;
			test.differential = 0;
			test.test_status = NOTEST;
		}
		else
		{
            //test.test_stat = 0;
			test.p_value = p_val;
			test.value_1 = 0;
			test.value_2 = 0;
			test.differential = js;
			test.test_status = enough_reads ? OK : NOTEST;
        }
	}
	else // we won't even bother with the JS-based testing in LOWDATA cases.
	{
        if (prev_status == NUMERIC_OK && curr_status == NUMERIC_OK && 
            prev_abundance.has_member_with_status(NUMERIC_LOW_DATA) == false &&
            curr_abundance.has_member_with_status(NUMERIC_LOW_DATA) == false)
            test.test_status = NOTEST;
        else if (prev_status == NUMERIC_FAIL || curr_status == NUMERIC_FAIL)
            test.test_status = FAIL;
        else
            test.test_status = LOWDATA;
            
		test.test_stat = 0;
		test.p_value = 0.0;
		test.differential = 0.0;
	}
    
    return test;
}

// This performs within-group tests on a set of isoforms or a set of TSS groups.
// This is a way of looking for meaningful differential parental expression between the same isoform/tss/cds/gene, differential parental splicing, or differential parental differential
// promoter use.
SampleDifference get_ds_tests(const AlleleAbundanceGroup& prev_abundance,
                              const AlleleAbundanceGroup& curr_abundance,
//                              SampleDiffs& diff_tests,
                              bool enough_reads,
							  bool prev_parent,
							  bool curr_parent,
							  bool doTest=true)
{	
	//the assignment in prev_parent,curr_parent: false(0) = paternal;true (1) = maternal
	const string& name = curr_abundance.description();
	
	SampleDifference test;
    
	test.test_status = NOTEST;
	
	AbundanceStatus prev_status,curr_status;
	double prev_num_fragments,curr_num_fragments;
	vector<ublas::vector<double> > sample_kappas;
	ublas::vector<double> curr_kappas(curr_abundance.abundances().size());
	ublas::vector<double> prev_kappas(prev_abundance.abundances().size());
	bool prev_member_with_status,curr_member_with_status;
	if(!doTest){
		test.p_value = 1;
        test.value_1 = 0;
        test.value_2 = 0;
        test.differential = 0;
        test.test_status = NOALLELETEST;
		return test;
	}
	
	if(!curr_parent)
	{
		curr_status = curr_abundance.paternal_status();
		curr_num_fragments = curr_abundance.num_paternal_fragments();
		for (size_t i = 0; i < curr_abundance.abundances().size(); ++i)
		{
			curr_kappas(i) = curr_abundance.abundances()[i]->paternal_kappa();
		}
		curr_member_with_status = curr_abundance.has_paternal_member_with_status(NUMERIC_LOW_DATA);
	}
	else
	{
		curr_status = curr_abundance.maternal_status();
		curr_num_fragments = curr_abundance.num_maternal_fragments();
		for (size_t i = 0; i < curr_abundance.abundances().size(); ++i)
		{
			curr_kappas(i) = curr_abundance.abundances()[i]->maternal_kappa();
		}
		curr_member_with_status = curr_abundance.has_maternal_member_with_status(NUMERIC_LOW_DATA);
	}
	if(!prev_parent)
	{
		prev_status = prev_abundance.paternal_status();
		prev_num_fragments = prev_abundance.num_paternal_fragments();
		for (size_t i = 0; i < prev_abundance.abundances().size(); ++i)
		{
			prev_kappas(i) = prev_abundance.abundances()[i]->paternal_kappa();
		}
		prev_member_with_status = prev_abundance.has_paternal_member_with_status(NUMERIC_LOW_DATA);
	}
	else
	{
		prev_status = prev_abundance.maternal_status();
		prev_num_fragments = prev_abundance.num_maternal_fragments();
		for (size_t i = 0; i < prev_abundance.abundances().size(); ++i)
		{
			prev_kappas(i) = prev_abundance.abundances()[i]->maternal_kappa();
		}
		prev_member_with_status = prev_abundance.has_maternal_member_with_status(NUMERIC_LOW_DATA);
	}	
    
    //if (prev_abundance.abundances().size() == 1 ||
    //    (prev_status == NUMERIC_OK && prev_num_fragments == 0) ||
    //    (curr_status == NUMERIC_OK && curr_num_fragments == 0))
	if ((prev_status == NUMERIC_OK && prev_num_fragments == 0) ||
        (curr_status == NUMERIC_OK && curr_num_fragments == 0))
    {
        test.p_value = 1;
        test.value_1 = 0;
        test.value_2 = 0;
        test.differential = 0;
        test.test_status = NOTEST;
    }
	//else if (prev_abundance.abundances().size() > 1 &&
        /*prev_abundance.has_member_with_status(NUMERIC_LOW_DATA) == false &&
        filtered_curr.has_member_with_status(NUMERIC_LOW_DATA) == false &&*/ 
	//  prev_status == NUMERIC_OK && prev_num_fragments > 0 &&
	//curr_status == NUMERIC_OK && curr_num_fragments > 0)
	else if (prev_status == NUMERIC_OK && prev_num_fragments > 0 &&
			 curr_status == NUMERIC_OK && curr_num_fragments > 0)
        /*prev_abundance.has_member_with_status(NUMERIC_LOW_DATA) == false &&
        filtered_curr.has_member_with_status(NUMERIC_LOW_DATA) == false &&*/ 
	{
		sample_kappas.push_back(prev_kappas);
		sample_kappas.push_back(curr_kappas);
		
        double js = 0.0;
        double p_val = 1.0;
        
        bool success;
        success = test_js(prev_abundance, curr_abundance, js, p_val, prev_parent, curr_parent);
        
		if (js == 0.0 || success == false)
		{
			test.test_stat = 0;
			test.p_value = 1.0;
			test.value_1 = 0;
			test.value_2 = 0;
			test.differential = 0;
			test.test_status = NOTEST;
		}
		else
		{
            //test.test_stat = 0;
			test.p_value = p_val;
			test.value_1 = 0;
			test.value_2 = 0;
			test.differential = js;
			test.test_status = enough_reads ? OK : NOTEST;
        }
	}
	else // we won't even bother with the JS-based testing in LOWDATA cases.
	{
        if (prev_status == NUMERIC_OK && curr_status == NUMERIC_OK && 
            prev_member_with_status == false &&
            curr_member_with_status == false)
            test.test_status = NOTEST;
        else if (prev_status == NUMERIC_FAIL || curr_status == NUMERIC_FAIL)
            test.test_status = FAIL;
        else
            test.test_status = LOWDATA;
            
		test.test_stat = 0;
		test.p_value = 0.0;
		test.differential = 0.0;
	}
    
    return test;
}

string make_ref_tag(const string& ref, char classcode)
{
	char tag_buf[1024];

	sprintf(tag_buf, 
			"%s(%c)",
			ref.c_str(),
			classcode);
	
	return string(tag_buf);
}

string bundle_locus_tag(const RefSequenceTable& rt, 
						const HitBundle& bundle)
{
	char locus_buf[1024];
	RefID bundle_chr_id = bundle.ref_id();
	assert (bundle_chr_id != 0);
	const char* chr_name = rt.get_name(bundle_chr_id);
	
	sprintf(locus_buf, 
			"%s:%d-%d",
			chr_name,
			bundle.left(),
			bundle.right());
	
	return string(locus_buf);
}

void sample_abundance_worker(const string& locus_tag,
                             const set<shared_ptr<ReadGroupProperties const> >& rg_props,
                             SampleAbundances& sample,
                             HitBundle* sample_bundle,
                             bool perform_cds_analysis,
                             bool perform_tss_analysis)
{
    vector<shared_ptr<Abundance> > abundances;
    
    foreach(shared_ptr<Scaffold> s, sample_bundle->ref_scaffolds())
    {
        TranscriptAbundance* pT = new TranscriptAbundance;
        pT->transfrag(s);
        shared_ptr<Abundance> ab(pT);
        ab->description(s->annotated_trans_id());
        ab->locus_tag(locus_tag);
        abundances.push_back(ab);
    }

    sample.transcripts = AbundanceGroup(abundances);
    sample.transcripts.init_rg_props(rg_props);
    
    vector<MateHit> hits_in_cluster;
    
    if (sample_bundle->hits().size() < (size_t)max_frags_per_bundle)
    {
		get_alignments_from_scaffolds(sample.transcripts.abundances(),
                                      hits_in_cluster);
    
		// Compute the individual transcript FPKMs via each sample's 
        // AbundanceGroup for this locus.
        
        sample.transcripts.calculate_abundance(hits_in_cluster);
	}
    else
    {
		foreach(shared_ptr<Abundance>  ab, abundances)
        {
            ab->status(NUMERIC_HI_DATA);

            CountPerReplicateTable cpr;
            FPKMPerReplicateTable fpr;
            StatusPerReplicateTable spr;
            for (set<shared_ptr<ReadGroupProperties const> >::const_iterator itr = rg_props.begin();
                 itr != rg_props.end();
                 ++itr)
            {
                cpr[*itr] = 0;
                fpr[*itr] = 0;
                spr[*itr] = NUMERIC_HI_DATA;
            }
            ab->num_fragments_by_replicate(cpr);
            ab->FPKM_by_replicate(fpr);
            ab->status_by_replicate(spr);
        }
        
    }
	
    // Cluster transcripts by gene_id
    vector<AbundanceGroup> transcripts_by_gene_id;
    cluster_transcripts<ConnectByAnnotatedGeneId>(sample.transcripts,
                                                  transcripts_by_gene_id);
    
	foreach(AbundanceGroup& ab_group, transcripts_by_gene_id)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
	
    sample.genes = transcripts_by_gene_id;
    
    if (perform_cds_analysis)
    {
        // Cluster transcripts by CDS
        vector<AbundanceGroup> transcripts_by_cds;
        ublas::matrix<double> cds_gamma_cov;
        ublas::matrix<double> cds_count_cov;
        ublas::matrix<double> cds_iterated_exp_count_cov;
        ublas::matrix<double> cds_fpkm_cov;
        vector<Eigen::VectorXd> cds_assigned_counts;
        
        vector<bool> mask(sample.transcripts.abundances().size(), true);
        for (size_t i = 0; i < sample.transcripts.abundances().size(); ++i)
        {
            if (*(sample.transcripts.abundances()[i]->protein_id().begin()) == "")
            {
                mask[i] = false;
            }
        }
        
        AbundanceGroup trans_with_p_id; 
        sample.transcripts.filter_group(mask, trans_with_p_id);
        
        cluster_transcripts<ConnectByAnnotatedProteinId>(trans_with_p_id,
                                                         transcripts_by_cds,
                                                         &cds_gamma_cov,
                                                         &cds_iterated_exp_count_cov,
                                                         &cds_count_cov,
                                                         &cds_fpkm_cov,
                                                         &cds_assigned_counts);
        
        foreach(AbundanceGroup& ab_group, transcripts_by_cds)
        {
            ab_group.locus_tag(locus_tag);
            set<string> protein_ids = ab_group.protein_id();
            assert (protein_ids.size() == 1);
            string desc = *(protein_ids.begin()); 
            //if (desc != "")
            //{
            assert (desc != "");
            ab_group.description(*(protein_ids.begin()));
            //}
        }
        
        sample.cds = transcripts_by_cds;
        
        // Group the CDS clusters by gene
        vector<shared_ptr<Abundance> > cds_abundances;
        double max_cds_mass_variance = 0.0; 
        set<shared_ptr<ReadGroupProperties const> > rg_props;
        foreach (AbundanceGroup& ab_group, sample.cds)
        {
            //if (ab_group.description() != "")
            {
                cds_abundances.push_back(shared_ptr<Abundance>(new AbundanceGroup(ab_group)));
                max_cds_mass_variance = max(ab_group.max_mass_variance(), max_cds_mass_variance);
                rg_props.insert(ab_group.rg_props().begin(), ab_group.rg_props().end()); 
            }
        }
        AbundanceGroup cds(cds_abundances,
                           cds_gamma_cov,
                           cds_iterated_exp_count_cov,
                           cds_count_cov,
                           cds_fpkm_cov,
                           max_cds_mass_variance,
                           rg_props,
                           cds_assigned_counts);
        
        vector<AbundanceGroup> cds_by_gene;
        
        cluster_transcripts<ConnectByAnnotatedGeneId>(cds,
                                                      cds_by_gene);
        
        foreach(AbundanceGroup& ab_group, cds_by_gene)
        {
            ab_group.locus_tag(locus_tag);
            set<string> gene_ids = ab_group.gene_id();
            assert (gene_ids.size() == 1);
            ab_group.description(*(gene_ids.begin()));
        }
        
        sample.gene_cds = cds_by_gene;
    }
    
    if (perform_tss_analysis)
    {
        // Cluster transcripts by start site (TSS)
        vector<AbundanceGroup> transcripts_by_tss;
        
        ublas::matrix<double> tss_gamma_cov;
        ublas::matrix<double> tss_count_cov;
        ublas::matrix<double> tss_iterated_exp_count_cov;
        ublas::matrix<double> tss_fpkm_cov;
        vector<Eigen::VectorXd> tss_assigned_counts;
        
        vector<bool> mask(sample.transcripts.abundances().size(), true);
        for (size_t i = 0; i < sample.transcripts.abundances().size(); ++i)
        {
            if (*(sample.transcripts.abundances()[i]->tss_id().begin()) == "")
            {
                mask[i] = false;
            }
        }
        
        AbundanceGroup trans_with_tss; 
        sample.transcripts.filter_group(mask, trans_with_tss);
        
        cluster_transcripts<ConnectByAnnotatedTssId>(trans_with_tss,
                                                     transcripts_by_tss,
                                                     &tss_gamma_cov,
                                                     &tss_iterated_exp_count_cov,
                                                     &tss_count_cov,
                                                     &tss_fpkm_cov,
                                                     &tss_assigned_counts);
        
        
        foreach(AbundanceGroup& ab_group, transcripts_by_tss)
        {
            ab_group.locus_tag(locus_tag);
            set<string> tss_ids = ab_group.tss_id();
            assert (tss_ids.size() == 1);
            string desc = *(tss_ids.begin()); 
            assert (desc != "");
            ab_group.description(*(tss_ids.begin()));
            
        }
        
        sample.primary_transcripts = transcripts_by_tss;
        double max_tss_mass_variance = 0.0;
        
        // Group TSS clusters by gene
        vector<shared_ptr<Abundance> > primary_transcript_abundances;
        set<shared_ptr<ReadGroupProperties const> > rg_props;
        foreach (AbundanceGroup& ab_group, sample.primary_transcripts)
        {
            primary_transcript_abundances.push_back(shared_ptr<Abundance>(new AbundanceGroup(ab_group)));
            max_tss_mass_variance = max(max_tss_mass_variance, ab_group.max_mass_variance());
            rg_props.insert(ab_group.rg_props().begin(), ab_group.rg_props().end());
        }
        
        AbundanceGroup primary_transcripts(primary_transcript_abundances,
                                           tss_gamma_cov,
                                           tss_iterated_exp_count_cov,
                                           tss_count_cov,
                                           tss_fpkm_cov,
                                           max_tss_mass_variance,
                                           rg_props,
                                           tss_assigned_counts);
        
        vector<AbundanceGroup> primary_transcripts_by_gene;
        
        cluster_transcripts<ConnectByAnnotatedGeneId>(primary_transcripts,
                                                      primary_transcripts_by_gene);
        
        foreach(AbundanceGroup& ab_group, primary_transcripts_by_gene)
        {
            ab_group.locus_tag(locus_tag);
            set<string> gene_ids = ab_group.gene_id();
//            if (gene_ids.size() > 1)
//            {
//                foreach (string st, gene_ids)
//                {
//                    fprintf(stderr, "%s\n", st.c_str());
//                }
//                ab_group.gene_id();
//            }
            assert (gene_ids.size() == 1);
            ab_group.description(*(gene_ids.begin()));
        }
        
        sample.gene_primary_transcripts = primary_transcripts_by_gene;
    }
}

void sample_abundance_worker(const string& locus_tag,
                             const set<shared_ptr<ReadGroupProperties const> >& rg_props,
                             SampleAlleleAbundances& sample,
                             HitBundle* sample_bundle,
                             bool perform_cds_analysis,
                             bool perform_tss_analysis)
{
    vector<shared_ptr<Abundance> > abundances;
    
    foreach(shared_ptr<Scaffold> s, sample_bundle->ref_scaffolds())
    {
        AlleleTranscriptAbundance* pT = new AlleleTranscriptAbundance;
        pT->transfrag(s);
		pT->set_allele_informative();
        shared_ptr<Abundance> ab(pT);
        ab->description(s->annotated_trans_id());
        ab->locus_tag(locus_tag);
        abundances.push_back(ab);
    }
	sample.transcripts = AlleleAbundanceGroup(abundances);
    sample.transcripts.init_rg_props(rg_props);
    
    vector<MateHit> hits_in_cluster;
    
    if (sample_bundle->hits().size() < (size_t)max_frags_per_bundle)
    {
		get_alignments_from_scaffolds(sample.transcripts.abundances(),
                                      hits_in_cluster);
        
        // Compute the individual transcript FPKMs via each sample's 
        // AbundanceGroup for this locus.
        
        sample.transcripts.calculate_abundance(hits_in_cluster);
    }
    else
    {
		foreach(shared_ptr<Abundance>  ab, abundances)
        {
            ab->paternal_status(NUMERIC_HI_DATA);
			ab->maternal_status(NUMERIC_HI_DATA);
            CountPerReplicateTable paternal_cpr,maternal_cpr;
            FPKMPerReplicateTable paternal_fpr,maternal_fpr;
            StatusPerReplicateTable paternal_spr,maternal_spr;
            for (set<shared_ptr<ReadGroupProperties const> >::const_iterator itr = rg_props.begin();
                 itr != rg_props.end();
                 ++itr)
            {
                paternal_cpr[*itr] = 0;
				maternal_cpr[*itr] = 0;
                paternal_fpr[*itr] = 0;
				maternal_fpr[*itr] = 0;
                paternal_spr[*itr] = NUMERIC_HI_DATA;
				maternal_spr[*itr] = NUMERIC_HI_DATA;
            }
            ab->num_paternal_fragments_by_replicate(paternal_cpr);
			ab->num_maternal_fragments_by_replicate(maternal_cpr);
            ab->paternal_FPKM_by_replicate(paternal_fpr);
			ab->maternal_FPKM_by_replicate(maternal_fpr);
            ab->paternal_status_by_replicate(paternal_spr);
			ab->maternal_status_by_replicate(maternal_spr);
        }
        
    }
    
    // Cluster transcripts by gene_id
    vector<AlleleAbundanceGroup> transcripts_by_gene_id;
    cluster_transcripts<ConnectByAnnotatedGeneId>(sample.transcripts,
                                                  transcripts_by_gene_id);
    
	foreach(AlleleAbundanceGroup& ab_group, transcripts_by_gene_id)
    {
        ab_group.locus_tag(locus_tag);
        set<string> gene_ids = ab_group.gene_id();
        assert (gene_ids.size() == 1);
        ab_group.description(*(gene_ids.begin()));
    }
	
    sample.genes = transcripts_by_gene_id;
    
    if (perform_cds_analysis)
    {
        // Cluster transcripts by CDS
        vector<AlleleAbundanceGroup> transcripts_by_cds;
        ublas::matrix<double> cds_gamma_cov;
        ublas::matrix<double> cds_count_cov;
        ublas::matrix<double> cds_iterated_exp_count_cov;
        ublas::matrix<double> cds_fpkm_cov;
        vector<Eigen::VectorXd> paternal_cds_assigned_counts;
		vector<Eigen::VectorXd> maternal_cds_assigned_counts;
        
        vector<bool> mask(sample.transcripts.abundances().size(), true);
        for (size_t i = 0; i < sample.transcripts.abundances().size(); ++i)
        {
            if (*(sample.transcripts.abundances()[i]->protein_id().begin()) == "")
            {
                mask[i] = false;
            }
        }
        
        AlleleAbundanceGroup trans_with_p_id; 
        sample.transcripts.filter_group(mask, trans_with_p_id);
        
        cluster_transcripts<ConnectByAnnotatedProteinId>(trans_with_p_id,
                                                         transcripts_by_cds,
                                                         &cds_gamma_cov,
                                                         &cds_iterated_exp_count_cov,
                                                         &cds_count_cov,
                                                         &cds_fpkm_cov,
                                                         &paternal_cds_assigned_counts,
														 &maternal_cds_assigned_counts);
        
        foreach(AlleleAbundanceGroup& ab_group, transcripts_by_cds)
        {
            ab_group.locus_tag(locus_tag);
            set<string> protein_ids = ab_group.protein_id();
            assert (protein_ids.size() == 1);
            string desc = *(protein_ids.begin()); 
            //if (desc != "")
            //{
            assert (desc != "");
            ab_group.description(*(protein_ids.begin()));
            //}
        }
        
        sample.cds = transcripts_by_cds;
        
        // Group the CDS clusters by gene
        vector<shared_ptr<Abundance> > cds_abundances;
        double max_cds_mass_variance = 0.0; 
        set<shared_ptr<ReadGroupProperties const> > rg_props;
        foreach (AlleleAbundanceGroup& ab_group, sample.cds)
        {
            //if (ab_group.description() != "")
            {
                cds_abundances.push_back(shared_ptr<Abundance>(new AlleleAbundanceGroup(ab_group)));
                max_cds_mass_variance = max(ab_group.max_mass_variance(), max_cds_mass_variance);
                rg_props.insert(ab_group.rg_props().begin(), ab_group.rg_props().end()); 
            }
        }
        AlleleAbundanceGroup cds(cds_abundances,
								 cds_gamma_cov,
								 cds_iterated_exp_count_cov,
								 cds_count_cov,
								 cds_fpkm_cov,
								 max_cds_mass_variance,
								 rg_props,
								 paternal_cds_assigned_counts,
								 maternal_cds_assigned_counts);
        
        vector<AlleleAbundanceGroup> cds_by_gene;
        
        cluster_transcripts<ConnectByAnnotatedGeneId>(cds,
                                                      cds_by_gene);
        
        foreach(AlleleAbundanceGroup& ab_group, cds_by_gene)
        {
            ab_group.locus_tag(locus_tag);
            set<string> gene_ids = ab_group.gene_id();
            assert (gene_ids.size() == 1);
            ab_group.description(*(gene_ids.begin()));
        }
        
        sample.gene_cds = cds_by_gene;
    }
	if (perform_tss_analysis)
    {
        // Cluster transcripts by start site (TSS)
        vector<AlleleAbundanceGroup> transcripts_by_tss;
        
        ublas::matrix<double> tss_gamma_cov;
        ublas::matrix<double> tss_count_cov;
        ublas::matrix<double> tss_iterated_exp_count_cov;
        ublas::matrix<double> tss_fpkm_cov;
        vector<Eigen::VectorXd> paternal_tss_assigned_counts;
		vector<Eigen::VectorXd> maternal_tss_assigned_counts;
        
        vector<bool> mask(sample.transcripts.abundances().size(), true);
        for (size_t i = 0; i < sample.transcripts.abundances().size(); ++i)
        {
            if (*(sample.transcripts.abundances()[i]->tss_id().begin()) == "")
            {
                mask[i] = false;
            }
        }
        
        AlleleAbundanceGroup trans_with_tss; 
        sample.transcripts.filter_group(mask, trans_with_tss);
        
        cluster_transcripts<ConnectByAnnotatedTssId>(trans_with_tss,
                                                     transcripts_by_tss,
                                                     &tss_gamma_cov,
                                                     &tss_iterated_exp_count_cov,
                                                     &tss_count_cov,
                                                     &tss_fpkm_cov,
                                                     &paternal_tss_assigned_counts,
													 &maternal_tss_assigned_counts);
        
        
        foreach(AlleleAbundanceGroup& ab_group, transcripts_by_tss)
        {
            ab_group.locus_tag(locus_tag);
            set<string> tss_ids = ab_group.tss_id();
            assert (tss_ids.size() == 1);
            string desc = *(tss_ids.begin()); 
            assert (desc != "");
            ab_group.description(*(tss_ids.begin()));
            
        }
        
        sample.primary_transcripts = transcripts_by_tss;
        double max_tss_mass_variance = 0.0;
        
        // Group TSS clusters by gene
        vector<shared_ptr<Abundance> > primary_transcript_abundances;
        set<shared_ptr<ReadGroupProperties const> > rg_props;
        foreach (AlleleAbundanceGroup& ab_group, sample.primary_transcripts)
        {
            primary_transcript_abundances.push_back(shared_ptr<Abundance>(new AlleleAbundanceGroup(ab_group)));
            max_tss_mass_variance = max(max_tss_mass_variance, ab_group.max_mass_variance());
            rg_props.insert(ab_group.rg_props().begin(), ab_group.rg_props().end());
        }
        
        AlleleAbundanceGroup primary_transcripts(primary_transcript_abundances,
												 tss_gamma_cov,
												 tss_iterated_exp_count_cov,
												 tss_count_cov,
												 tss_fpkm_cov,
												 max_tss_mass_variance,
												 rg_props,
												 paternal_tss_assigned_counts,
												 maternal_tss_assigned_counts);
        
        vector<AlleleAbundanceGroup> primary_transcripts_by_gene;
        
        cluster_transcripts<ConnectByAnnotatedGeneId>(primary_transcripts,
                                                      primary_transcripts_by_gene);
        
        foreach(AlleleAbundanceGroup& ab_group, primary_transcripts_by_gene)
        {
            ab_group.locus_tag(locus_tag);
            set<string> gene_ids = ab_group.gene_id();
//            if (gene_ids.size() > 1)
//            {
//                foreach (string st, gene_ids)
//                {
//                    fprintf(stderr, "%s\n", st.c_str());
//                }
//                ab_group.gene_id();
//            }
            assert (gene_ids.size() == 1);
            ab_group.description(*(gene_ids.begin()));
        }
        
        sample.gene_primary_transcripts = primary_transcripts_by_gene;
    }
}

struct LocusVarianceInfo
{
    int factory_id;
    double count_mean;
    double count_empir_var;
    double locus_count_fitted_var;
    double isoform_fitted_var_sum;
    double cross_replicate_js;
    int num_transcripts;
    double bayes_gamma_trace;
    double empir_gamma_trace;
    vector<double> gamma;
    vector<double> gamma_var;
    vector<double> gamma_bootstrap_var;
    vector<string> transcript_ids;
    vector<double> count_sharing;
    double locus_salient_frags;
    double locus_total_frags;

};

#if ENABLE_THREADS
mutex variance_info_lock; // don't modify the above struct without locking here
#endif

vector<LocusVarianceInfo> locus_variance_info_table;

// We'll use this tracking table to collect per replicate counts for each 
// transcript, so we can re-fit the variance model.
FPKMTrackingTable transcript_count_tracking; 

void sample_worker(const RefSequenceTable& rt,
                   ReplicatedBundleFactory& sample_factory,
                   shared_ptr<SampleAbundances> abundance,
                   size_t factory_id,
                   shared_ptr<TestLauncher> launcher)
{
#if ENABLE_THREADS
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif
    
    HitBundle bundle;
    bool non_empty = sample_factory.next_bundle(bundle);
    
    if (!non_empty || (bias_run && bundle.ref_scaffolds().size() != 1)) // Only learn on single isoforms
    {
#if !ENABLE_THREADS
        // If Cuffdiff was built without threads, we need to manually invoke 
        // the testing functor, which will check to see if all the workers
        // are done, and if so, perform the cross sample testing.
        launcher->abundance_avail(locus_tag, abundance, factory_id);
        launcher->test_finished_loci();
        //launcher();
#endif
    	return;
    }

    abundance->cluster_mass = bundle.mass();
    
    char bundle_label_buf[2048];
    sprintf(bundle_label_buf, 
            "%s:%d-%d", 
            rt.get_name(bundle.ref_id()),
            bundle.left(),
            bundle.right());
    string locus_tag = bundle_label_buf;
    
    launcher->register_locus(locus_tag);
    
    abundance->locus_tag = locus_tag;
    
    bool perform_cds_analysis = false;
    bool perform_tss_analysis = false;

    foreach(shared_ptr<Scaffold> s, bundle.ref_scaffolds())
    {
        if (s->annotated_tss_id() != "")
        {
            perform_tss_analysis = final_est_run;
        }
        if (s->annotated_protein_id() != "")
        {
            perform_cds_analysis = final_est_run;
        }
    }

    set<shared_ptr<ReadGroupProperties const> > rg_props;
    for (size_t i = 0; i < sample_factory.factories().size(); ++i)
    {
        shared_ptr<BundleFactory> bf = sample_factory.factories()[i];
        rg_props.insert(bf->read_group_properties());
    }
    
    sample_abundance_worker(boost::cref(locus_tag),
                            boost::cref(rg_props),
                            boost::ref(*abundance),
                            &bundle,
                            perform_cds_analysis,
                            perform_tss_analysis);
    
#if ENABLE_THREADS
    variance_info_lock.lock();
#endif
    
    //fprintf(stderr, "\nTesting in %s (%d total tests)\n", locus_tag.c_str(), total_tests);
    
	// Add all the transcripts, CDS groups, TSS groups, and genes to their
    // respective FPKM tracking table.  Whether this is a time series or an
    // all pairs comparison, we should be calculating and reporting FPKMs for 
    // all objects in all samples
    
//    ///////////////////////////////////////////////
//    shared_ptr<MassDispersionModel const> disperser = sample_factory.mass_dispersion_model();
//    pair<double, double> locus_mv = disperser->get_raw_mean_and_var(locus_tag);
//    if (locus_mv.first != 0 && locus_mv.second != 0)
//    {
//        LocusVarianceInfo info;
//        info.factory_id = factory_id;
//        info.count_mean = locus_mv.first;
//        info.count_empir_var = locus_mv.second;
//        info.locus_count_fitted_var = disperser->scale_mass_variance(info.count_mean);
//        
//        double total_iso_scaled_var = 0.0;
//        
//        const AbundanceGroup& ab_group = abundance->transcripts;
//        info.locus_total_frags = ab_group.total_frags();
//        info.locus_salient_frags = ab_group.salient_frags();
//        //double group_counts = ab_group.total_frags();
//        ublas::matrix<double> cov = ab_group.iterated_count_cov();
//        if (ab_group.num_fragments())
//            cov /= ab_group.num_fragments();
//        
//        double total_length = 0.0;
//        for (unsigned i = 0; i < ab_group.abundances().size(); ++ i) 
//        {
//            total_length += ab_group.abundances()[i]->effective_length();
//        }
//        
////        if (total_length)
////        {
////            for (unsigned i = 0; i < ab_group.abundances().size(); ++ i) 
////            {
////                fprintf(stderr, 
////                        "%lg, %lg, %lg\n", 
////                        _abundances[i]->gamma(), 
////                        _abundances[i]->effective_length()/total_length, 
////                        log2(_abundances[i]->gamma()/(_abundances[i]->effective_length()/total_length)));
////            }
////        }
//        
//		for (size_t i = 0; i < ab_group.abundances().size(); ++i)
//		{
//            
////            double count_var = cov(i,i);
////            double max_count_covar = 0.0;
////            size_t max_covar_idx = 0.0;
////            for (size_t j = 0; j < cov.size1(); ++j)
////            {
////                if (j != i && abs(cov(i,j)) > max_count_covar)
////                {
////                    max_count_covar = abs(cov(i,j));
////                    max_covar_idx = j;
////                }
////            }
//            double count_sharing = 0.0;
////            if (cov(i,i) != 0 && cov(max_covar_idx,max_covar_idx) != 0)
////                count_sharing = -1.0 * cov(i,max_covar_idx) / sqrt(cov(i,i) * cov(max_covar_idx,max_covar_idx));
//            
//            
//            if (total_length)
//                count_sharing = log2(ab_group.abundances()[i]->gamma()/(ab_group.abundances()[i]->effective_length()/total_length));
//            
//            shared_ptr<Abundance> ab = ab_group.abundances()[i];
//            double scaled_var = disperser->scale_mass_variance(ab->num_fragments());
//			total_iso_scaled_var += scaled_var;
//            info.gamma.push_back(ab->gamma());
//            info.gamma_var.push_back(ab_group.gamma_cov()(i,i));
//            info.count_sharing.push_back(count_sharing);
//            info.transcript_ids.push_back(ab->description());
//		}
//
//        
//        const ublas::matrix<double>& gamma_cov = ab_group.gamma_cov();
//        info.bayes_gamma_trace = 0;
//        info.empir_gamma_trace = 0;
//        for (size_t i = 0; i < ab_group.abundances().size(); ++i)
//        {
//            //for (size_t j = 0; j < ab_group.abundances().size(); ++j)
//            {
//                info.bayes_gamma_trace += gamma_cov(i,i);
//            }
//        }
//
//        
//        info.cross_replicate_js = 30;
//        //assert (abundance->cluster_mass == locus_mv.first);
//        //assert (total_iso_scaled_var >= info.count_mean);
//        
//        info.isoform_fitted_var_sum = total_iso_scaled_var;
//        info.num_transcripts = ab_group.abundances().size();
////        info.bayes_gamma_trace = 0;
////        info.empir_gamma_trace = 0;
//        locus_variance_info_table.push_back(info);
//    }
    
#if ENABLE_THREADS
    variance_info_lock.unlock();
#endif
    
    ///////////////////////////////////////////////
    
    
    foreach(shared_ptr<Scaffold> ref_scaff,  bundle.ref_scaffolds())
    {
        ref_scaff->clear_hits();
    }
    
    launcher->abundance_avail(locus_tag, abundance, factory_id);
    launcher->test_finished_loci();
    
#if !ENABLE_THREADS
    // If Cuffdiff was built without threads, we need to manually invoke 
    // the testing functor, which will check to see if all the workers
    // are done, and if so, perform the cross sample testing.
    //launcher->test_finished_loci();
#endif
}

//nimrod
void allele_sample_worker(const RefSequenceTable& rt,
                   ReplicatedBundleFactory& sample_factory,
                   shared_ptr<SampleAlleleAbundances> abundance,
                   size_t factory_id,
                   shared_ptr<AlleleTestLauncher> launcher)
{
#if ENABLE_THREADS
	boost::this_thread::at_thread_exit(decr_pool_count);
#endif
    
    HitBundle bundle;
    bool non_empty = sample_factory.next_bundle(bundle);
    
    if (!non_empty || (bias_run && bundle.ref_scaffolds().size() != 1)) // Only learn on single isoforms
    {
#if !ENABLE_THREADS
        // If Cuffdiff was built without threads, we need to manually invoke 
        // the testing functor, which will check to see if all the workers
        // are done, and if so, perform the cross sample testing.
        launcher->abundance_avail(locus_tag, abundance, factory_id);
        launcher->test_finished_loci();
        //launcher();
#endif
    	return;
    }
    
    abundance->cluster_mass = bundle.mass();
    
    char bundle_label_buf[2048];
    sprintf(bundle_label_buf, 
            "%s:%d-%d", 
            rt.get_name(bundle.ref_id()),
            bundle.left(),
            bundle.right());
    string locus_tag = bundle_label_buf;
    
    launcher->register_locus(locus_tag);
    
    abundance->locus_tag = locus_tag;
    
    bool perform_cds_analysis = false;
    bool perform_tss_analysis = false;

    foreach(shared_ptr<Scaffold> s, bundle.ref_scaffolds())
    {
        if (s->annotated_tss_id() != "")
        {
            perform_tss_analysis = final_est_run;
        }
        if (s->annotated_protein_id() != "")
        {
            perform_cds_analysis = final_est_run;
        }
    }

    set<shared_ptr<ReadGroupProperties const> > rg_props;
    for (size_t i = 0; i < sample_factory.factories().size(); ++i)
    {
        shared_ptr<BundleFactory> bf = sample_factory.factories()[i];
        rg_props.insert(bf->read_group_properties());
    }
    
    sample_abundance_worker(boost::cref(locus_tag),
                            boost::cref(rg_props),
                            boost::ref(*abundance),
                            &bundle,
                            perform_cds_analysis,
                            perform_tss_analysis);
	
#if ENABLE_THREADS
    variance_info_lock.lock();
#endif
    
    //fprintf(stderr, "\nTesting in %s (%d total tests)\n", locus_tag.c_str(), total_tests);
    
	// Add all the transcripts, CDS groups, TSS groups, and genes to their
    // respective FPKM tracking table.  Whether this is a time series or an
    // all pairs comparison, we should be calculating and reporting FPKMs for 
    // all objects in all samples
    
//    ///////////////////////////////////////////////
//    shared_ptr<MassDispersionModel const> disperser = sample_factory.mass_dispersion_model();
//    pair<double, double> locus_mv = disperser->get_raw_mean_and_var(locus_tag);
//    if (locus_mv.first != 0 && locus_mv.second != 0)
//    {
//        LocusVarianceInfo info;
//        info.factory_id = factory_id;
//        info.count_mean = locus_mv.first;
//        info.count_empir_var = locus_mv.second;
//        info.locus_count_fitted_var = disperser->scale_mass_variance(info.count_mean);
//        
//        double total_iso_scaled_var = 0.0;
//        
//        const AbundanceGroup& ab_group = abundance->transcripts;
//        info.locus_total_frags = ab_group.total_frags();
//        info.locus_salient_frags = ab_group.salient_frags();
//        //double group_counts = ab_group.total_frags();
//        ublas::matrix<double> cov = ab_group.iterated_count_cov();
//        if (ab_group.num_fragments())
//            cov /= ab_group.num_fragments();
//        
//        double total_length = 0.0;
//        for (unsigned i = 0; i < ab_group.abundances().size(); ++ i) 
//        {
//            total_length += ab_group.abundances()[i]->effective_length();
//        }
//        
////        if (total_length)
////        {
////            for (unsigned i = 0; i < ab_group.abundances().size(); ++ i) 
////            {
////                fprintf(stderr, 
////                        "%lg, %lg, %lg\n", 
////                        _abundances[i]->gamma(), 
////                        _abundances[i]->effective_length()/total_length, 
////                        log2(_abundances[i]->gamma()/(_abundances[i]->effective_length()/total_length)));
////            }
////        }
//        
//		for (size_t i = 0; i < ab_group.abundances().size(); ++i)
//		{
//            
////            double count_var = cov(i,i);
////            double max_count_covar = 0.0;
////            size_t max_covar_idx = 0.0;
////            for (size_t j = 0; j < cov.size1(); ++j)
////            {
////                if (j != i && abs(cov(i,j)) > max_count_covar)
////                {
////                    max_count_covar = abs(cov(i,j));
////                    max_covar_idx = j;
////                }
////            }
//            double count_sharing = 0.0;
////            if (cov(i,i) != 0 && cov(max_covar_idx,max_covar_idx) != 0)
////                count_sharing = -1.0 * cov(i,max_covar_idx) / sqrt(cov(i,i) * cov(max_covar_idx,max_covar_idx));
//            
//            
//            if (total_length)
//                count_sharing = log2(ab_group.abundances()[i]->gamma()/(ab_group.abundances()[i]->effective_length()/total_length));
//            
//            shared_ptr<Abundance> ab = ab_group.abundances()[i];
//            double scaled_var = disperser->scale_mass_variance(ab->num_fragments());
//			total_iso_scaled_var += scaled_var;
//            info.gamma.push_back(ab->gamma());
//            info.gamma_var.push_back(ab_group.gamma_cov()(i,i));
//            info.count_sharing.push_back(count_sharing);
//            info.transcript_ids.push_back(ab->description());
//		}
//
//        
//        const ublas::matrix<double>& gamma_cov = ab_group.gamma_cov();
//        info.bayes_gamma_trace = 0;
//        info.empir_gamma_trace = 0;
//        for (size_t i = 0; i < ab_group.abundances().size(); ++i)
//        {
//            //for (size_t j = 0; j < ab_group.abundances().size(); ++j)
//            {
//                info.bayes_gamma_trace += gamma_cov(i,i);
//            }
//        }
//
//        
//        info.cross_replicate_js = 30;
//        //assert (abundance->cluster_mass == locus_mv.first);
//        //assert (total_iso_scaled_var >= info.count_mean);
//        
//        info.isoform_fitted_var_sum = total_iso_scaled_var;
//        info.num_transcripts = ab_group.abundances().size();
////        info.bayes_gamma_trace = 0;
////        info.empir_gamma_trace = 0;
//        locus_variance_info_table.push_back(info);
//    }
    
#if ENABLE_THREADS
    variance_info_lock.unlock();
#endif
    
    ///////////////////////////////////////////////
    
    
    foreach(shared_ptr<Scaffold> ref_scaff,  bundle.ref_scaffolds())
    {
        ref_scaff->clear_hits();
    }
    
    launcher->abundance_avail(locus_tag, abundance, factory_id);
	launcher->test_finished_loci();
    
#if !ENABLE_THREADS
    // If Cuffdiff was built without threads, we need to manually invoke 
    // the testing functor, which will check to see if all the workers
    // are done, and if so, perform the cross sample testing.
    //launcher->test_finished_loci();
#endif
}

void dump_locus_variance_info(const string& filename)
{
#if ENABLE_THREADS
    variance_info_lock.lock();
#endif
    
    FILE* fdump = fopen(filename.c_str(), "w");
    
    fprintf(fdump, 
            "condition\tdescription\tlocus_counts\tempir_var\tlocus_fit_var\tsum_iso_fit_var\tcross_replicate_js\tnum_transcripts\tbayes_gamma_trace\tempir_gamma_trace\tcount_mean\tgamma_var\tlocus_salient_frags\tlocus_total_frags\tcount_sharing\n");
    foreach (LocusVarianceInfo& L, locus_variance_info_table)
    {
        for (size_t i = 0; i < L.gamma.size(); ++i)
        {
            fprintf(fdump, "%d\t%s\t%lg\t%lg\t%lg\t%lg\t%lg\t%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n", L.factory_id, L.transcript_ids[i].c_str(), L.count_mean, L.count_empir_var, L.locus_count_fitted_var, L.isoform_fitted_var_sum, L.cross_replicate_js, L.num_transcripts, L.bayes_gamma_trace, L.empir_gamma_trace,L.gamma[i],L.gamma_var[i], L.locus_salient_frags, L.locus_total_frags, L.count_sharing[i]);
        }
        
    }
    
#if ENABLE_THREADS
    variance_info_lock.unlock();
#endif
}

void filter_group_for_js_testing(vector<vector<AbundanceGroup> >& source_groups)
{
    if (source_groups.empty())
        return;
    
    // iterate over transcript groups
    for (size_t ab_group_idx = 0; ab_group_idx < source_groups[0].size(); ++ab_group_idx)
    {
        vector<bool> to_keep(source_groups[0][ab_group_idx].abundances().size(), false);
        // iterate over samples
        for (size_t sample_idx = 0;  sample_idx < source_groups.size(); ++sample_idx)
        {
            // iterate over each member of the current abundance group
            AbundanceGroup& ab_group = source_groups[sample_idx][ab_group_idx];
            for (size_t ab_idx = 0; ab_idx < ab_group.abundances().size(); ++ab_idx)
            {
                const Abundance& ab = *(ab_group.abundances()[ab_idx]);
                if (ab.num_fragments() && ab.effective_length())
                {
                    double frags_per_kb = ab.num_fragments() / (ab.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        to_keep[ab_idx] = true;
                }
            }
        }
        
        // Now that we know which ones we want to keep, get rid of the rest
        for (size_t sample_idx = 0;  sample_idx < source_groups.size(); ++sample_idx)
        {
            AbundanceGroup& ab_group = source_groups[sample_idx][ab_group_idx];
            AbundanceGroup f = ab_group;
            ab_group.filter_group(to_keep, f);
            ab_group = f;
        }
    }
}

void filter_group_for_js_testing(vector<vector<AlleleAbundanceGroup> >& source_groups)
{
    if (source_groups.empty())
        return;
    
    // iterate over transcript groups
    for (size_t ab_group_idx = 0; ab_group_idx < source_groups[0].size(); ++ab_group_idx)
    {
        vector<bool> to_keep(source_groups[0][ab_group_idx].abundances().size(), false);
        // iterate over samples
        for (size_t sample_idx = 0;  sample_idx < source_groups.size(); ++sample_idx)
        {
            // iterate over each member of the current abundance group
            AlleleAbundanceGroup& ab_group = source_groups[sample_idx][ab_group_idx];
            for (size_t ab_idx = 0; ab_idx < ab_group.abundances().size(); ++ab_idx)
            {
                const Abundance& ab = *(ab_group.abundances()[ab_idx]);
				if (ab.is_allele_informative() && (ab.num_paternal_fragments()+ab.num_maternal_fragments()) && (ab.paternal_effective_length() || ab.maternal_effective_length()))
                {
                    double paternal_frags_per_kb = ab.num_paternal_fragments() / (ab.paternal_effective_length() / 1000.0);
					double maternal_frags_per_kb = ab.num_maternal_fragments() / (ab.maternal_effective_length() / 1000.0);
					if (paternal_frags_per_kb >= min_read_count || maternal_frags_per_kb >= min_read_count)
                        to_keep[ab_idx] = true;
                }
            }
        }
        
        // Now that we know which ones we want to keep, get rid of the rest
        for (size_t sample_idx = 0;  sample_idx < source_groups.size(); ++sample_idx)
        {
            AlleleAbundanceGroup& ab_group = source_groups[sample_idx][ab_group_idx];
            AlleleAbundanceGroup f = ab_group;
            ab_group.filter_group(to_keep, f);
            ab_group = f;
        }
    }
}

bool group_has_record_above_thresh(const AbundanceGroup& ab_group)
{
    for (size_t ab_idx = 0; ab_idx < ab_group.abundances().size(); ++ab_idx)
    {
        const Abundance& ab = *(ab_group.abundances()[ab_idx]);
        if (ab.num_fragments() && ab.effective_length())
        {
            double frags_per_kb = ab.num_fragments() / (ab.effective_length() / 1000.0);
            if (frags_per_kb >= min_read_count)
                return true;
        }
    }
    return false;
}

bool group_has_record_above_thresh_paternal(const AlleleAbundanceGroup& ab_group)
{
    for (size_t ab_idx = 0; ab_idx < ab_group.abundances().size(); ++ab_idx)
    {
        const Abundance& ab = *(ab_group.abundances()[ab_idx]);
		if (ab.num_paternal_fragments() && ab.paternal_effective_length())
        {
            double frags_per_kb = ab.num_paternal_fragments() / (ab.paternal_effective_length() / 1000.0);
			if (frags_per_kb >= min_read_count)
                return true;
        }
    }
	return false;
}

bool group_has_record_above_thresh_maternal(const AlleleAbundanceGroup& ab_group)
{
    for (size_t ab_idx = 0; ab_idx < ab_group.abundances().size(); ++ab_idx)
    {
        const Abundance& ab = *(ab_group.abundances()[ab_idx]);
        if (ab.num_maternal_fragments() && ab.maternal_effective_length())
        {
            double frags_per_kb = ab.num_maternal_fragments() / (ab.maternal_effective_length() / 1000.0);
            if (frags_per_kb >= min_read_count)
                return true;
        }
    }
    return false;
}


bool is_badly_fit(const Abundance& ab)
{
    double pooled_fpkm = ab.FPKM();
    double pooled_fpkm_var = ab.FPKM_variance();
    if (pooled_fpkm == 0 || pooled_fpkm_var <= 0)
        return false;
    normal norm(pooled_fpkm, sqrt(pooled_fpkm_var));
    FPKMPerReplicateTable fpkm_by_rep = ab.FPKM_by_replicate();
    for (FPKMPerReplicateTable::const_iterator f_itr = fpkm_by_rep.begin();
         f_itr != fpkm_by_rep.end(); ++f_itr)
    {
        double rep_fpkm = f_itr->second;
        double p_value = cdf(norm, rep_fpkm);
        if (p_value < min_outlier_p)
            return true;
    }
    return false;
}

bool is_badly_fit_paternal(const Abundance& ab)
{
    double pooled_fpkm = ab.paternal_FPKM();
    double pooled_fpkm_var = ab.paternal_FPKM_variance();
    if (pooled_fpkm == 0 || pooled_fpkm_var <= 0)
        return false;
    normal norm(pooled_fpkm, sqrt(pooled_fpkm_var));
    FPKMPerReplicateTable fpkm_by_rep = ab.paternal_FPKM_by_replicate();
	int rep = 0; //DEBUG
	for (FPKMPerReplicateTable::const_iterator f_itr = fpkm_by_rep.begin();
         f_itr != fpkm_by_rep.end(); ++f_itr)
    {
		rep += 1;//DEBUG
		double rep_fpkm = f_itr->second;
        double p_value = cdf(norm, rep_fpkm);
        if (p_value < min_outlier_p/2.0)
			return true;
	}
    return false;
}

bool is_badly_fit_maternal(const Abundance& ab)
{
    double pooled_fpkm = ab.maternal_FPKM();
    double pooled_fpkm_var = ab.maternal_FPKM_variance();
    if (pooled_fpkm == 0 || pooled_fpkm_var <= 0)
        return false;
    normal norm(pooled_fpkm, sqrt(pooled_fpkm_var));
    FPKMPerReplicateTable fpkm_by_rep = ab.maternal_FPKM_by_replicate();
	int rep = 0; //DEBUG
    for (FPKMPerReplicateTable::const_iterator f_itr = fpkm_by_rep.begin();
         f_itr != fpkm_by_rep.end(); ++f_itr)
    {
		rep += 1;//DEBUG
        double rep_fpkm = f_itr->second;
        double p_value = cdf(norm, rep_fpkm);
        if (p_value < min_outlier_p/2.0)
            return true;
	}
    return false;
}

bool group_has_record_badly_fit(const AbundanceGroup& ab_group)
{
    for (size_t ab_idx = 0; ab_idx < ab_group.abundances().size(); ++ab_idx)
    {
        const Abundance& ab = *(ab_group.abundances()[ab_idx]);
        if (ab.num_fragments() && ab.effective_length())
        {
            if (is_badly_fit(ab))
                return true;
        }
    }
    return false;
}

bool group_has_paternal_record_badly_fit(const AlleleAbundanceGroup& ab_group)
{
    for (size_t ab_idx = 0; ab_idx < ab_group.abundances().size(); ++ab_idx)
    {
        const Abundance& ab = *(ab_group.abundances()[ab_idx]);
		if (ab.num_paternal_fragments() && ab.paternal_effective_length())
        {
			if (is_badly_fit_paternal(ab))
                return true;
        }
    }
    return false;
}

bool group_has_maternal_record_badly_fit(const AlleleAbundanceGroup& ab_group)
{
    for (size_t ab_idx = 0; ab_idx < ab_group.abundances().size(); ++ab_idx)
    {
        const Abundance& ab = *(ab_group.abundances()[ab_idx]);
        if (ab.num_maternal_fragments() && ab.maternal_effective_length())
        {
            if (is_badly_fit_maternal(ab))
                return true;
        }
    }
    return false;
}

int total_tests = 0;
void test_differential(const string& locus_tag,
					   const vector<shared_ptr<SampleAbundances> >& samples,
					   Tests& tests,
					   Tracking& tracking,
                       bool samples_are_time_series)
{
	if (samples.empty())
		return;
    
    if (no_differential == true)
    {
//#if ENABLE_THREADS
//        test_storage_lock.unlock();
//#endif
        return;
    }
	
    vector<vector<AbundanceGroup> > filtered_primary_trans_groups;
    vector<vector<AbundanceGroup> > filtered_promoter_groups;
    vector<vector<AbundanceGroup> > filtered_cds_groups;

    for (size_t i = 0; i < samples.size(); ++i)
    {
        filtered_primary_trans_groups.push_back(samples[i]->primary_transcripts);
        filtered_promoter_groups.push_back(samples[i]->gene_primary_transcripts);
        filtered_cds_groups.push_back(samples[i]->gene_cds);
    }
    
    filter_group_for_js_testing(filtered_primary_trans_groups);
    filter_group_for_js_testing(filtered_promoter_groups);
    filter_group_for_js_testing(filtered_cds_groups);
    
    
    // Perform pairwise significance testing between samples. If this is a
    // time series, only test between successive pairs of samples, as supplied 
    // by the user.
	for (size_t i = 1; i < samples.size(); ++i)
	{
		//bool multi_transcript_locus = samples[i]->transcripts.abundances().size() > 1;
		
        int sample_to_start_test_against = 0;
        if (samples_are_time_series)
            sample_to_start_test_against = i - 1;
        
        for (size_t j = sample_to_start_test_against; j < i; ++j)
        {
//            bool enough_reads = (samples[i]->cluster_mass >= min_read_count ||
//                                 samples[j]->cluster_mass >= min_read_count);
            assert (samples[i]->transcripts.abundances().size() == 
                    samples[j]->transcripts.abundances().size());
            for (size_t k = 0; k < samples[i]->transcripts.abundances().size(); ++k)
            {
                const Abundance& curr_abundance = *(samples[j]->transcripts.abundances()[k]);
                const Abundance& prev_abundance = *(samples[i]->transcripts.abundances()[k]);
                const string& desc = curr_abundance.description();
                FPKMTrackingTable::iterator itr = tracking.isoform_fpkm_tracking.find(desc);
                assert (itr != tracking.isoform_fpkm_tracking.end());
                
                bool enough_reads = false;
                if (curr_abundance.num_fragments() && curr_abundance.effective_length())
                {
                    double frags_per_kb = curr_abundance.num_fragments() / (curr_abundance.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        enough_reads = true;
                }
                if (prev_abundance.num_fragments() && prev_abundance.effective_length())
                {
                    double frags_per_kb = prev_abundance.num_fragments() / (prev_abundance.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        enough_reads = true;
                }
                
                if (enough_reads)
                {
                    if (is_badly_fit(curr_abundance) || is_badly_fit(prev_abundance))
                        enough_reads = false;
                }

                SampleDifference test;
				test = get_de_tests(desc,
                                      itr->second.fpkm_series[j], 
                                      itr->second.fpkm_series[i],
                                      //tests.isoform_de_tests[i][j],
                                      enough_reads);
#if ENABLE_THREADS
                test_storage_lock.lock();
#endif
                pair<SampleDiffs::iterator, bool> inserted; 
                inserted = tests.isoform_de_tests[i][j].insert(make_pair(desc,
                                                                     test)); 
                
                shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
                
                meta_data->gene_ids = curr_abundance.gene_id();
                meta_data->gene_names = curr_abundance.gene_name();
                meta_data->protein_ids = curr_abundance.protein_id();
                meta_data->locus_desc = curr_abundance.locus_tag();
                meta_data->description = curr_abundance.description();
                inserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
                test_storage_lock.unlock();
#endif
            }
            
            for (size_t k = 0; k < samples[i]->cds.size(); ++k)
            {
                const Abundance& curr_abundance = samples[i]->cds[k];
                const Abundance& prev_abundance = samples[j]->cds[k];
                
                const string& desc = curr_abundance.description();
                FPKMTrackingTable::iterator itr = tracking.cds_fpkm_tracking.find(desc);
                assert (itr != tracking.cds_fpkm_tracking.end());
                
                bool enough_reads = false;
                if (curr_abundance.num_fragments() && curr_abundance.effective_length())
                {
                    double frags_per_kb = curr_abundance.num_fragments() / (curr_abundance.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        enough_reads = true;
                }
                if (prev_abundance.num_fragments() && prev_abundance.effective_length())
                {
                    double frags_per_kb = prev_abundance.num_fragments() / (prev_abundance.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        enough_reads = true;
                }
                
                if (enough_reads)
                {
                    if (is_badly_fit(curr_abundance) || is_badly_fit(prev_abundance))
                        enough_reads = false;
                }

                
                SampleDifference test;
                test = get_de_tests(desc,
                                    itr->second.fpkm_series[j], 
                                    itr->second.fpkm_series[i],
                                    //tests.cds_de_tests[i][j],
                                    enough_reads);
#if ENABLE_THREADS
                test_storage_lock.lock();
#endif
                
                pair<SampleDiffs::iterator, bool> inserted; 
                inserted = tests.cds_de_tests[i][j].insert(make_pair(desc,
                                                                           test)); 
                
                shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
                
                meta_data->gene_ids = curr_abundance.gene_id();
                meta_data->gene_names = curr_abundance.gene_name();
                meta_data->protein_ids = curr_abundance.protein_id();
                meta_data->locus_desc = curr_abundance.locus_tag();
                meta_data->description = curr_abundance.description();
                inserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
                test_storage_lock.unlock();
#endif
            }
            
            for (size_t k = 0; k < samples[i]->primary_transcripts.size(); ++k)
            {
                const Abundance& curr_abundance = samples[i]->primary_transcripts[k];
                const Abundance& prev_abundance = samples[j]->primary_transcripts[k];
                
                const string& desc = curr_abundance.description();
                FPKMTrackingTable::iterator itr = tracking.tss_group_fpkm_tracking.find(desc);
                assert (itr != tracking.tss_group_fpkm_tracking.end());
                
                bool enough_reads = false;
                if (curr_abundance.num_fragments() && curr_abundance.effective_length())
                {
                    double frags_per_kb = curr_abundance.num_fragments() / (curr_abundance.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        enough_reads = true;
                }
                if (prev_abundance.num_fragments() && prev_abundance.effective_length())
                {
                    double frags_per_kb = prev_abundance.num_fragments() / (prev_abundance.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        enough_reads = true;
                }
                
                if (enough_reads)
                {
                    if (is_badly_fit(curr_abundance) || is_badly_fit(prev_abundance))
                        enough_reads = false;
                }

                
                SampleDifference test;
                test = get_de_tests(desc,
                             itr->second.fpkm_series[j], 
                             itr->second.fpkm_series[i],
                             //tests.tss_group_de_tests[i][j],
                             enough_reads);
#if ENABLE_THREADS
                test_storage_lock.lock();
#endif
                pair<SampleDiffs::iterator, bool> inserted; 
                inserted = tests.tss_group_de_tests[i][j].insert(make_pair(desc,
                                                                      test)); 
                
                
                shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
                
                meta_data->gene_ids = curr_abundance.gene_id();
                meta_data->gene_names = curr_abundance.gene_name();
                meta_data->protein_ids = curr_abundance.protein_id();
                meta_data->locus_desc = curr_abundance.locus_tag();
                meta_data->description = curr_abundance.description();
                inserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
                test_storage_lock.unlock();
#endif
            }
            
            for (size_t k = 0; k < samples[i]->genes.size(); ++k)
            {
                const AbundanceGroup& curr_abundance = samples[i]->genes[k];
                const AbundanceGroup& prev_abundance = samples[j]->genes[k];
                const string& desc = curr_abundance.description();
                FPKMTrackingTable::iterator itr = tracking.gene_fpkm_tracking.find(desc);
                assert (itr != tracking.gene_fpkm_tracking.end());
                
                bool enough_reads = false;
                if (curr_abundance.num_fragments() && curr_abundance.effective_length())
                {
                    double frags_per_kb = curr_abundance.num_fragments() / (curr_abundance.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        enough_reads = true;
                }
                if (prev_abundance.num_fragments() && prev_abundance.effective_length())
                {
                    double frags_per_kb = prev_abundance.num_fragments() / (prev_abundance.effective_length() / 1000.0);
                    if (frags_per_kb >= min_read_count)
                        enough_reads = true;
                }
                    
                SampleDifference test;
                test = get_de_tests(desc,
                             itr->second.fpkm_series[j], 
                             itr->second.fpkm_series[i],
                             //tests.gene_de_tests[i][j],
                             enough_reads);
#if ENABLE_THREADS
                test_storage_lock.lock();
#endif
                pair<SampleDiffs::iterator, bool> inserted; 
                inserted = tests.gene_de_tests[i][j].insert(make_pair(desc,
                                                                      test)); 
                
                shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
                
                meta_data->gene_ids = curr_abundance.gene_id();
                meta_data->gene_names = curr_abundance.gene_name();
                meta_data->protein_ids = curr_abundance.protein_id();
                meta_data->locus_desc = curr_abundance.locus_tag();
                meta_data->description = curr_abundance.description();
                inserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
                test_storage_lock.unlock();
#endif
            }
            
            // Skip all the JS based testing for genes with an isoform switch?
            if (no_js_tests)
                continue;
            
            // FIXME: the code below will not properly test for differential
            // splicing/promoter use when a gene (e.g.) occupies two
            // disjoint bundles.  We need to store the covariance matrices (etc)
            // in the FPKMContexts to handle that case properly.
            
            // Differential promoter use
            for (size_t k = 0; k < samples[i]->gene_primary_transcripts.size(); ++k)
            {
                const AbundanceGroup& curr_abundance = filtered_promoter_groups[i][k];
                const AbundanceGroup& prev_abundance = filtered_promoter_groups[j][k];
                const string& desc = curr_abundance.description();
                
                bool enough_reads = (group_has_record_above_thresh(curr_abundance) &&
                                     group_has_record_above_thresh(prev_abundance) &&
                                     group_has_record_badly_fit(curr_abundance) == false &&
                                     group_has_record_badly_fit(prev_abundance) == false &&
                                     curr_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test && 
                                     prev_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test);
                SampleDifference test;
                test = get_ds_tests(prev_abundance, 
                                    curr_abundance,
                                    enough_reads);
                
                // The filtered group might be empty, so let's grab metadata from
                // the unfiltered group
                shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
                
                meta_data->gene_ids = samples[i]->gene_primary_transcripts[k].gene_id();
                meta_data->gene_names = samples[i]->gene_primary_transcripts[k].gene_name();
                meta_data->protein_ids = samples[i]->gene_primary_transcripts[k].protein_id();
                meta_data->locus_desc = samples[i]->gene_primary_transcripts[k].locus_tag();
                meta_data->description = samples[i]->gene_primary_transcripts[k].description();
                
                test.meta_data = meta_data;
                
#if ENABLE_THREADS
                test_storage_lock.lock();
#endif
                
                pair<SampleDiffs::iterator, bool> inserted;
                inserted = tests.diff_promoter_tests[i][j].insert(make_pair(desc,test)); 
                inserted.first->second = test;
                
#if ENABLE_THREADS
                test_storage_lock.unlock();
#endif
            }
            
            // Differential coding sequence output
            for (size_t k = 0; k < samples[i]->gene_cds.size(); ++k)
            {
                const AbundanceGroup& curr_abundance = filtered_cds_groups[i][k];
                const AbundanceGroup& prev_abundance = filtered_cds_groups[j][k];
                const string& desc = curr_abundance.description();
                
                bool enough_reads = (group_has_record_above_thresh(curr_abundance) &&
                                     group_has_record_above_thresh(prev_abundance) &&
                                     group_has_record_badly_fit(curr_abundance) == false &&
                                     group_has_record_badly_fit(prev_abundance) == false &&
                                     curr_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test && 
                                     prev_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test);                
                SampleDifference test;
                test = get_ds_tests(prev_abundance, 
                                    curr_abundance,
                                    enough_reads);
                
                // The filtered group might be empty, so let's grab metadata from
                // the unfiltered group
                shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
                
                meta_data->gene_ids = samples[i]->gene_cds[k].gene_id();
                meta_data->gene_names = samples[i]->gene_cds[k].gene_name();
                meta_data->protein_ids = samples[i]->gene_cds[k].protein_id();
                meta_data->locus_desc = samples[i]->gene_cds[k].locus_tag();
                meta_data->description = samples[i]->gene_cds[k].description();
                
                test.meta_data = meta_data;
                
#if ENABLE_THREADS
                test_storage_lock.lock();
#endif
                pair<SampleDiffs::iterator, bool> inserted;
                inserted = tests.diff_cds_tests[i][j].insert(make_pair(desc,test)); 
                inserted.first->second = test;
                
#if ENABLE_THREADS
                test_storage_lock.unlock();
#endif
            }
            
            // Differential splicing of primary transcripts
            for (size_t k = 0; k < samples[i]->primary_transcripts.size(); ++k)
            {
                const AbundanceGroup& curr_abundance = filtered_primary_trans_groups[i][k];
                const AbundanceGroup& prev_abundance = filtered_primary_trans_groups[j][k];
                const string& desc = curr_abundance.description();
                
                bool enough_reads = (group_has_record_above_thresh(curr_abundance) &&
                                     group_has_record_above_thresh(prev_abundance) &&
                                     group_has_record_badly_fit(curr_abundance) == false &&
                                     group_has_record_badly_fit(prev_abundance) == false &&
                                     curr_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test && 
                                     prev_abundance.FPKM_by_replicate().size() >= min_reps_for_js_test);
                
                SampleDifference test;
                test = get_ds_tests(prev_abundance, 
                                    curr_abundance,
                                    enough_reads);
                
                // The filtered group might be empty, so let's grab metadata from
                // the unfiltered group
                shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
                
                meta_data->gene_ids = samples[i]->primary_transcripts[k].gene_id();
                meta_data->gene_names = samples[i]->primary_transcripts[k].gene_name();
                meta_data->protein_ids = samples[i]->primary_transcripts[k].protein_id();
                meta_data->locus_desc = samples[i]->primary_transcripts[k].locus_tag();
                meta_data->description = samples[i]->primary_transcripts[k].description();
                
                test.meta_data = meta_data;
                
#if ENABLE_THREADS
                test_storage_lock.lock();
#endif
                pair<SampleDiffs::iterator, bool> inserted;
                inserted = tests.diff_splicing_tests[i][j].insert(make_pair(desc,test)); 
                inserted.first->second = test;
                
#if ENABLE_THREADS
                test_storage_lock.unlock();
#endif
            }
        }
	}
}

void test_differential(const string& locus_tag,
					   const vector<shared_ptr<SampleAlleleAbundances> >& samples,
					   AlleleTests& tests,
					   Tracking& tracking,
					   bool samples_are_time_series)
{
	if (samples.empty())
		return;
    
    if (no_differential == true)
    {
//#if ENABLE_THREADS
//        test_storage_lock.unlock();
//#endif
        return;
    }
	
    vector<vector<AlleleAbundanceGroup> > filtered_primary_trans_groups;
    vector<vector<AlleleAbundanceGroup> > filtered_promoter_groups;
    vector<vector<AlleleAbundanceGroup> > filtered_cds_groups;

    for (size_t i = 0; i < samples.size(); ++i)
    {
        filtered_primary_trans_groups.push_back(samples[i]->primary_transcripts);
        filtered_promoter_groups.push_back(samples[i]->gene_primary_transcripts);
        filtered_cds_groups.push_back(samples[i]->gene_cds);
    }
	filter_group_for_js_testing(filtered_primary_trans_groups);
    filter_group_for_js_testing(filtered_promoter_groups);
    filter_group_for_js_testing(filtered_cds_groups);
    
    // Perform pairwise significance testing between samples. If this is a
    // time series, only test between successive pairs of samples, as supplied 
    // by the user.
	if(samples.size() > 1)
	{//expand all cuffdiff tests to pat_vs_pat, pat_vs_mat, mat_vs_mat, mat_vs_pat between any 2 samples
		for (size_t i = 1; i < samples.size(); ++i)
		{
			//bool multi_transcript_locus = samples[i]->transcripts.abundances().size() > 1;
			
			int sample_to_start_test_against = 0;
			if (samples_are_time_series)
				sample_to_start_test_against = i - 1;
			
			for (size_t j = sample_to_start_test_against; j < i; ++j)
			{
//            bool enough_reads = (samples[i]->cluster_mass >= min_read_count ||
//                                 samples[j]->cluster_mass >= min_read_count);
				assert (samples[i]->transcripts.abundances().size() == 
						samples[j]->transcripts.abundances().size());
				for (size_t k = 0; k < samples[i]->transcripts.abundances().size(); ++k)
				{
					const Abundance& curr_abundance = *(samples[j]->transcripts.abundances()[k]);
					const Abundance& prev_abundance = *(samples[i]->transcripts.abundances()[k]);
					const string& desc = curr_abundance.description();
					FPKMTrackingTable::iterator itr = tracking.isoform_fpkm_tracking.find(desc);
					assert (itr != tracking.isoform_fpkm_tracking.end());
										
					bool curr_paternal_enough_reads = false;
					bool curr_maternal_enough_reads = false;
					bool prev_paternal_enough_reads = false;
					bool prev_maternal_enough_reads = false;
					
					if ((curr_abundance.num_paternal_fragments()+curr_abundance.num_maternal_fragments()) && (curr_abundance.paternal_effective_length() || curr_abundance.maternal_effective_length()))
					{
						double paternal_frags_per_kb = curr_abundance.num_paternal_fragments() / (curr_abundance.paternal_effective_length() / 1000.0);
						double maternal_frags_per_kb = curr_abundance.num_maternal_fragments() / (curr_abundance.maternal_effective_length() / 1000.0);
						if (paternal_frags_per_kb >= min_read_count)
							curr_paternal_enough_reads = true;
						if (maternal_frags_per_kb >= min_read_count)
							curr_maternal_enough_reads = true;
					}
					if ((prev_abundance.num_paternal_fragments()+prev_abundance.num_maternal_fragments()) && (prev_abundance.paternal_effective_length() || prev_abundance.maternal_effective_length()))
					{
						double paternal_frags_per_kb = prev_abundance.num_paternal_fragments() / (prev_abundance.paternal_effective_length() / 1000.0);
						double maternal_frags_per_kb = prev_abundance.num_maternal_fragments() / (prev_abundance.maternal_effective_length() / 1000.0);
						if (paternal_frags_per_kb >= min_read_count)
							prev_paternal_enough_reads = true;
						if (maternal_frags_per_kb >= min_read_count)
							prev_maternal_enough_reads = true;
					}
					
					if (curr_paternal_enough_reads)
					{
						if (is_badly_fit_paternal(curr_abundance))
							curr_paternal_enough_reads = false;
					}
					if (curr_maternal_enough_reads)
					{
						if (is_badly_fit_maternal(curr_abundance))
							curr_maternal_enough_reads = false;
					}
					if (prev_paternal_enough_reads)
					{
						if (is_badly_fit_paternal(prev_abundance))
							curr_paternal_enough_reads = false;
					}
					if (prev_maternal_enough_reads)
					{
						if (is_badly_fit_maternal(prev_abundance))
							prev_maternal_enough_reads = false;
					}
					
					SampleDifference paternalPaternalTest;
					paternalPaternalTest = get_de_tests("CurrentPaternal:PreviousPaternal;"+desc,
														itr->second.paternal_fpkm_series[j], 
														itr->second.paternal_fpkm_series[i],
														//tests.isoform_de_tests[i][j],
														(curr_paternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference paternalMaternalTest;
					paternalMaternalTest = get_de_tests("CurrentPaternal:PreviousMaternal;"+desc,
														itr->second.paternal_fpkm_series[j], 
														itr->second.maternal_fpkm_series[i],
														//tests.isoform_de_tests[i][j],
														(curr_paternal_enough_reads || prev_maternal_enough_reads));
					SampleDifference maternalPaternalTest;
					maternalPaternalTest = get_de_tests("CurrentMaternal:PreviousPaternal;"+desc,
														itr->second.maternal_fpkm_series[j], 
														itr->second.paternal_fpkm_series[i],
														//tests.isoform_de_tests[i][j],
														(curr_maternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalMaternalTest;
					maternalMaternalTest = get_de_tests("CurrentMaternal:PreviousMaternal;"+desc,
														itr->second.maternal_fpkm_series[j], 
														itr->second.maternal_fpkm_series[i],
														//tests.isoform_de_tests[i][j],
														(curr_maternal_enough_reads || prev_maternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
#if ENABLE_THREADS
					test_storage_lock.lock();
#endif
					pair<SampleDiffs::iterator, bool> paternalPaternalInserted; 
					paternalPaternalInserted = tests.paternal_paternal_isoform_de_tests[i][j].insert(make_pair("CurrentPaternal:PreviousPaternal;"+desc,
																							 paternalPaternalTest)); 
					pair<SampleDiffs::iterator, bool> paternalMaternalInserted; 
					paternalMaternalInserted = tests.paternal_maternal_isoform_de_tests[i][j].insert(make_pair("CurrentPaternal:PreviousMaternal;"+desc,
																							 paternalMaternalTest)); 
					pair<SampleDiffs::iterator, bool> maternalPaternalInserted; 
					maternalPaternalInserted = tests.maternal_paternal_isoform_de_tests[i][j].insert(make_pair("CurrentMaternal:PreviousPaternal;"+desc,
																							 maternalPaternalTest)); 
					pair<SampleDiffs::iterator, bool> maternalMaternalInserted; 
					maternalMaternalInserted = tests.maternal_maternal_isoform_de_tests[i][j].insert(make_pair("CurrentMaternal:PreviousMaternal;"+desc,
																							 maternalMaternalTest)); 
					
					shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
					
					meta_data->gene_ids = curr_abundance.gene_id();
					meta_data->gene_names = curr_abundance.gene_name();
					meta_data->protein_ids = curr_abundance.protein_id();
					meta_data->locus_desc = curr_abundance.locus_tag();
					meta_data->description = curr_abundance.description();
					paternalPaternalInserted.first->second.meta_data = meta_data;
					paternalMaternalInserted.first->second.meta_data = meta_data;
					maternalPaternalInserted.first->second.meta_data = meta_data;
					maternalMaternalInserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
					test_storage_lock.unlock();
#endif
				}
				
				for (size_t k = 0; k < samples[i]->cds.size(); ++k)
				{
					const Abundance& curr_abundance = samples[i]->cds[k];
					const Abundance& prev_abundance = samples[j]->cds[k];
					
					const string& desc = curr_abundance.description();
					FPKMTrackingTable::iterator itr = tracking.cds_fpkm_tracking.find(desc);
					assert (itr != tracking.cds_fpkm_tracking.end());
										
					bool curr_paternal_enough_reads = false;
					bool curr_maternal_enough_reads = false;
					bool prev_paternal_enough_reads = false;
					bool prev_maternal_enough_reads = false;
					
					if ((curr_abundance.num_paternal_fragments()+curr_abundance.num_maternal_fragments()) && (curr_abundance.paternal_effective_length() || curr_abundance.maternal_effective_length()))
					{
						double paternal_frags_per_kb = curr_abundance.num_paternal_fragments() / (curr_abundance.paternal_effective_length() / 1000.0);
						double maternal_frags_per_kb = curr_abundance.num_maternal_fragments() / (curr_abundance.maternal_effective_length() / 1000.0);
						if (paternal_frags_per_kb >= min_read_count)
							curr_paternal_enough_reads = true;
						if (maternal_frags_per_kb >= min_read_count)
							curr_maternal_enough_reads = true;
					}
					if ((prev_abundance.num_paternal_fragments()+prev_abundance.num_maternal_fragments()) && (prev_abundance.paternal_effective_length() || prev_abundance.maternal_effective_length()))
					{
						double paternal_frags_per_kb = prev_abundance.num_paternal_fragments() / (prev_abundance.paternal_effective_length() / 1000.0);
						double maternal_frags_per_kb = prev_abundance.num_maternal_fragments() / (prev_abundance.maternal_effective_length() / 1000.0);
						if (paternal_frags_per_kb >= min_read_count)
							prev_paternal_enough_reads = true;
						if (maternal_frags_per_kb >= min_read_count)
							prev_maternal_enough_reads = true;
					}
					
					if (curr_paternal_enough_reads)
					{
						if (is_badly_fit_paternal(curr_abundance))
							curr_paternal_enough_reads = false;
					}
					if (curr_maternal_enough_reads)
					{
						if (is_badly_fit_maternal(curr_abundance))
							curr_maternal_enough_reads = false;
					}
					if (prev_paternal_enough_reads)
					{
						if (is_badly_fit_paternal(prev_abundance))
							curr_paternal_enough_reads = false;
					}
					if (prev_maternal_enough_reads)
					{
						if (is_badly_fit_maternal(prev_abundance))
							prev_maternal_enough_reads = false;
					}
					
					SampleDifference paternalPaternalTest;
					paternalPaternalTest = get_de_tests("CurrentPaternal:PreviousPaternal;"+desc,
														itr->second.paternal_fpkm_series[j], 
														itr->second.paternal_fpkm_series[i],
														//tests.cds_de_tests[i][j],
														(curr_paternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference paternalMaternalTest;
					paternalMaternalTest = get_de_tests("CurrentPaternal:PreviousMaternal;"+desc,
														itr->second.paternal_fpkm_series[j], 
														itr->second.maternal_fpkm_series[i],
														//tests.cds_de_tests[i][j],
														(curr_paternal_enough_reads || prev_maternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalPaternalTest;
					paternalPaternalTest = get_de_tests("CurrentMaternal:PreviousPaternal;"+desc,
														itr->second.maternal_fpkm_series[j], 
														itr->second.paternal_fpkm_series[i],
														//tests.cds_de_tests[i][j],
														(curr_maternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalMaternalTest;
					paternalPaternalTest = get_de_tests("CurrentMaternal:PreviousMaternal;"+desc,
														itr->second.maternal_fpkm_series[j], 
														itr->second.maternal_fpkm_series[i],
														//tests.cds_de_tests[i][j],
														(curr_paternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
#if ENABLE_THREADS
					test_storage_lock.lock();
#endif
					
					pair<SampleDiffs::iterator, bool> paternalPaternalInserted; 
					paternalPaternalInserted = tests.paternal_paternal_cds_de_tests[i][j].insert(make_pair("CurrentPaternal:PreviousPaternal;"+desc,
																						 paternalPaternalTest)); 
					pair<SampleDiffs::iterator, bool> paternalMaternalInserted; 
					paternalMaternalInserted = tests.paternal_maternal_cds_de_tests[i][j].insert(make_pair("CurrentPaternal:PreviousMaternal;"+desc,
																						 paternalMaternalTest)); 
					pair<SampleDiffs::iterator, bool> maternalPaternalInserted; 
					maternalPaternalInserted = tests.maternal_paternal_cds_de_tests[i][j].insert(make_pair("CurrentMaternal:PreviousPaternal;"+desc,
																						 maternalPaternalTest)); 
					pair<SampleDiffs::iterator, bool> maternalMaternalInserted; 
					maternalMaternalInserted = tests.maternal_maternal_cds_de_tests[i][j].insert(make_pair("CurrentMaternal:PreviousMaternal;"+desc,
																						 maternalMaternalTest)); 
					
					shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
					
					meta_data->gene_ids = curr_abundance.gene_id();
					meta_data->gene_names = curr_abundance.gene_name();
					meta_data->protein_ids = curr_abundance.protein_id();
					meta_data->locus_desc = curr_abundance.locus_tag();
					meta_data->description = curr_abundance.description();
					paternalPaternalInserted.first->second.meta_data = meta_data;
					paternalMaternalInserted.first->second.meta_data = meta_data;
					maternalPaternalInserted.first->second.meta_data = meta_data;
					maternalMaternalInserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
					test_storage_lock.unlock();
#endif
				}
				
				for (size_t k = 0; k < samples[i]->primary_transcripts.size(); ++k)
				{
					const Abundance& curr_abundance = samples[i]->primary_transcripts[k];
					const Abundance& prev_abundance = samples[j]->primary_transcripts[k];
					
					const string& desc = curr_abundance.description();
					FPKMTrackingTable::iterator itr = tracking.tss_group_fpkm_tracking.find(desc);
					assert (itr != tracking.tss_group_fpkm_tracking.end());
					
					bool curr_paternal_enough_reads = false;
					bool curr_maternal_enough_reads = false;
					bool prev_paternal_enough_reads = false;
					bool prev_maternal_enough_reads = false;
					
					
					if ((curr_abundance.num_paternal_fragments()+curr_abundance.num_maternal_fragments()) && (curr_abundance.paternal_effective_length() || curr_abundance.maternal_effective_length()))
					{
						double paternal_frags_per_kb = curr_abundance.num_paternal_fragments() / (curr_abundance.paternal_effective_length() / 1000.0);
						double maternal_frags_per_kb = curr_abundance.num_maternal_fragments() / (curr_abundance.maternal_effective_length() / 1000.0);
						if (paternal_frags_per_kb >= min_read_count)
							curr_paternal_enough_reads = true;
						if (maternal_frags_per_kb >= min_read_count)
							curr_maternal_enough_reads = true;
					}
					if ((prev_abundance.num_paternal_fragments()+prev_abundance.num_maternal_fragments()) && (prev_abundance.paternal_effective_length() || prev_abundance.maternal_effective_length()))
					{
						double paternal_frags_per_kb = prev_abundance.num_paternal_fragments() / (prev_abundance.paternal_effective_length() / 1000.0);
						double maternal_frags_per_kb = prev_abundance.num_maternal_fragments() / (prev_abundance.maternal_effective_length() / 1000.0);
						if (paternal_frags_per_kb >= min_read_count)
							prev_paternal_enough_reads = true;
						if (maternal_frags_per_kb >= min_read_count)
							prev_maternal_enough_reads = true;
					}
					
					if (curr_paternal_enough_reads)
					{
						if (is_badly_fit_paternal(curr_abundance))
							curr_paternal_enough_reads = false;
					}
					if (curr_maternal_enough_reads)
					{
						if (is_badly_fit_maternal(curr_abundance))
							curr_maternal_enough_reads = false;
					}
					if (prev_paternal_enough_reads)
					{
						if (is_badly_fit_paternal(prev_abundance))
							curr_paternal_enough_reads = false;
					}
					if (prev_maternal_enough_reads)
					{
						if (is_badly_fit_maternal(prev_abundance))
							prev_maternal_enough_reads = false;
					}
					
					SampleDifference paternalPaternalTest;
					paternalPaternalTest = get_de_tests("CurrentPaternal:PreviousPaternal;"+desc,
														itr->second.paternal_fpkm_series[j], 
														itr->second.paternal_fpkm_series[i],
														//tests.tss_de_tests[i][j],
														(curr_paternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference paternalMaternalTest;
					paternalMaternalTest = get_de_tests("CurrentPaternal:PreviousMaternal;"+desc,
														itr->second.paternal_fpkm_series[j], 
														itr->second.maternal_fpkm_series[i],
														//tests.tss_de_tests[i][j],
														(curr_paternal_enough_reads || prev_maternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalPaternalTest;
					paternalPaternalTest = get_de_tests("CurrentMaternal:PreviousPaternal;"+desc,
														itr->second.maternal_fpkm_series[j], 
														itr->second.paternal_fpkm_series[i],
														//tests.tss_de_tests[i][j],
														(curr_maternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalMaternalTest;
					paternalPaternalTest = get_de_tests("CurrentMaternal:PreviousMaternal;"+desc,
														itr->second.maternal_fpkm_series[j], 
														itr->second.maternal_fpkm_series[i],
														//tests.tss_de_tests[i][j],
														(curr_paternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					
#if ENABLE_THREADS
					test_storage_lock.lock();
#endif
					pair<SampleDiffs::iterator, bool> paternalPaternalInserted; 
					paternalPaternalInserted = tests.paternal_paternal_tss_group_de_tests[i][j].insert(make_pair("CurrentPaternal:PreviousPaternal;"+desc,
																						 paternalPaternalTest)); 
					pair<SampleDiffs::iterator, bool> paternalMaternalInserted; 
					paternalMaternalInserted = tests.paternal_maternal_tss_group_de_tests[i][j].insert(make_pair("CurrentPaternal:PreviousMaternal;"+desc,
																						 paternalMaternalTest)); 
					pair<SampleDiffs::iterator, bool> maternalPaternalInserted; 
					maternalPaternalInserted = tests.maternal_paternal_tss_group_de_tests[i][j].insert(make_pair("CurrentMaternal:PreviousPaternal;"+desc,
																						 maternalPaternalTest)); 
					pair<SampleDiffs::iterator, bool> maternalMaternalInserted; 
					maternalMaternalInserted = tests.maternal_maternal_tss_group_de_tests[i][j].insert(make_pair("CurrentMaternal:PreviousMaternal;"+desc,
																						 maternalMaternalTest)); 
					
					shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
					
					meta_data->gene_ids = curr_abundance.gene_id();
					meta_data->gene_names = curr_abundance.gene_name();
					meta_data->protein_ids = curr_abundance.protein_id();
					meta_data->locus_desc = curr_abundance.locus_tag();
					meta_data->description = curr_abundance.description();
					paternalPaternalInserted.first->second.meta_data = meta_data;
					paternalMaternalInserted.first->second.meta_data = meta_data;
					maternalPaternalInserted.first->second.meta_data = meta_data;
					maternalMaternalInserted.first->second.meta_data = meta_data;
#if ENABLE_THREADS
					test_storage_lock.unlock();
#endif
				}
				
				for (size_t k = 0; k < samples[i]->genes.size(); ++k)
				{
					const AlleleAbundanceGroup& curr_abundance = samples[i]->genes[k];
					const AlleleAbundanceGroup& prev_abundance = samples[j]->genes[k];
					const string& desc = curr_abundance.description();
					FPKMTrackingTable::iterator itr = tracking.gene_fpkm_tracking.find(desc);
					assert (itr != tracking.gene_fpkm_tracking.end());
					
					bool curr_paternal_enough_reads = false;
					bool curr_maternal_enough_reads = false;
					bool prev_paternal_enough_reads = false;
					bool prev_maternal_enough_reads = false;
					
					if ((curr_abundance.num_paternal_fragments()+curr_abundance.num_maternal_fragments()) && (curr_abundance.paternal_effective_length() || curr_abundance.maternal_effective_length()))
					{
						double paternal_frags_per_kb = curr_abundance.num_paternal_fragments() / (curr_abundance.paternal_effective_length() / 1000.0);
						double maternal_frags_per_kb = curr_abundance.num_maternal_fragments() / (curr_abundance.maternal_effective_length() / 1000.0);
						if (paternal_frags_per_kb >= min_read_count)
							curr_paternal_enough_reads = true;
						if (maternal_frags_per_kb >= min_read_count)
							curr_maternal_enough_reads = true;
					}
					if ((prev_abundance.num_paternal_fragments()+prev_abundance.num_maternal_fragments()) && (prev_abundance.paternal_effective_length() || prev_abundance.maternal_effective_length()))
					{
						double paternal_frags_per_kb = prev_abundance.num_paternal_fragments() / (prev_abundance.paternal_effective_length() / 1000.0);
						double maternal_frags_per_kb = prev_abundance.num_maternal_fragments() / (prev_abundance.maternal_effective_length() / 1000.0);
						if (paternal_frags_per_kb >= min_read_count)
							prev_paternal_enough_reads = true;
						if (maternal_frags_per_kb >= min_read_count)
							prev_maternal_enough_reads = true;
					}
                    
					SampleDifference paternalPaternalTest;
					paternalPaternalTest = get_de_tests("CurrentPaternal:PreviousPaternal;"+desc,
														itr->second.paternal_fpkm_series[j], 
														itr->second.paternal_fpkm_series[i],
														//tests.gene_de_tests[i][j],
														(curr_paternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference paternalMaternalTest;
					paternalMaternalTest = get_de_tests("CurrentPaternal:PreviousMaternal;"+desc,
														itr->second.paternal_fpkm_series[j], 
														itr->second.maternal_fpkm_series[i],
														//tests.gene_de_tests[i][j],
														(curr_paternal_enough_reads || prev_maternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalPaternalTest;
					paternalPaternalTest = get_de_tests("CurrentMaternal:PreviousPaternal;"+desc,
														itr->second.maternal_fpkm_series[j], 
														itr->second.paternal_fpkm_series[i],
														//tests.gene_de_tests[i][j],
														(curr_maternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalMaternalTest;
					paternalPaternalTest = get_de_tests("CurrentMaternal:PreviousMaternal;"+desc,
														itr->second.maternal_fpkm_series[j], 
														itr->second.maternal_fpkm_series[i],
														//tests.gene_de_tests[i][j],
														(curr_paternal_enough_reads || prev_paternal_enough_reads),
														(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					
#if ENABLE_THREADS
					test_storage_lock.lock();
#endif
					pair<SampleDiffs::iterator, bool> paternalPaternalInserted; 
					paternalPaternalInserted = tests.paternal_paternal_gene_de_tests[i][j].insert(make_pair("CurrentPaternal:PreviousPaternal;"+desc,
																						  paternalPaternalTest)); 
					pair<SampleDiffs::iterator, bool> paternalMaternalInserted; 
					paternalMaternalInserted = tests.paternal_maternal_gene_de_tests[i][j].insert(make_pair("CurrentPaternal:PreviousMaternal;"+desc,
																						  paternalMaternalTest)); 
					pair<SampleDiffs::iterator, bool> maternalPaternalInserted; 
					maternalPaternalInserted = tests.maternal_paternal_gene_de_tests[i][j].insert(make_pair("CurrentMaternal:PreviousPaternal;"+desc,
																						  maternalPaternalTest)); 
					pair<SampleDiffs::iterator, bool> maternalMaternalInserted; 
					maternalMaternalInserted = tests.maternal_maternal_gene_de_tests[i][j].insert(make_pair("CurrentMaternal:PreviousMaternal;"+desc,
																						  maternalMaternalTest)); 
					
					shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
					
					meta_data->gene_ids = curr_abundance.gene_id();
					meta_data->gene_names = curr_abundance.gene_name();
					meta_data->protein_ids = curr_abundance.protein_id();
					meta_data->locus_desc = curr_abundance.locus_tag();
					meta_data->description = curr_abundance.description();
					paternalPaternalInserted.first->second.meta_data = meta_data;
					paternalMaternalInserted.first->second.meta_data = meta_data;
					maternalPaternalInserted.first->second.meta_data = meta_data;
					maternalMaternalInserted.first->second.meta_data = meta_data;
					
#if ENABLE_THREADS
					test_storage_lock.unlock();
#endif
				}
				
				// Skip all the JS based testing for genes with an isoform switch?
				if (no_js_tests)
					continue;
				
				// FIXME: the code below will not properly test for differential
				// splicing/promoter use when a gene (e.g.) occupies two
				// disjoint bundles.  We need to store the covariance matrices (etc)
				// in the FPKMContexts to handle that case properly.
				
				// Differential promoter use
				for (size_t k = 0; k < samples[i]->gene_primary_transcripts.size(); ++k)
				{
					const AlleleAbundanceGroup& curr_abundance = filtered_promoter_groups[i][k];
					const AlleleAbundanceGroup& prev_abundance = filtered_promoter_groups[j][k];
					const string& desc = curr_abundance.description();
					
					bool paternal_paternal_enough_reads = (group_has_record_above_thresh_paternal(curr_abundance) &&
														   group_has_record_above_thresh_paternal(prev_abundance) &&
														   group_has_paternal_record_badly_fit(curr_abundance) == false &&
														   group_has_paternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					bool paternal_maternal_enough_reads = (group_has_record_above_thresh_paternal(curr_abundance) &&
														   group_has_record_above_thresh_maternal(prev_abundance) &&
														   group_has_paternal_record_badly_fit(curr_abundance) == false &&
														   group_has_maternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					bool maternal_paternal_enough_reads = (group_has_record_above_thresh_maternal(curr_abundance) &&
														   group_has_record_above_thresh_paternal(prev_abundance) &&
														   group_has_maternal_record_badly_fit(curr_abundance) == false &&
														   group_has_paternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					bool maternal_maternal_enough_reads = (group_has_record_above_thresh_maternal(curr_abundance) &&
														   group_has_record_above_thresh_maternal(prev_abundance) &&
														   group_has_maternal_record_badly_fit(curr_abundance) == false &&
														   group_has_maternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					SampleDifference paternalPaternalTest;
					paternalPaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														paternal_paternal_enough_reads,
														false,false,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference paternalMaternalTest;
					paternalMaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														paternal_maternal_enough_reads,
														true,false,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalPaternalTest;
					maternalPaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														maternal_paternal_enough_reads,
														false,true,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalMaternalTest;
					maternalMaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														maternal_maternal_enough_reads,
														true,true,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					
					// The filtered group might be empty, so let's grab metadata from
					// the unfiltered group
					shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
					
					meta_data->gene_ids = samples[i]->gene_primary_transcripts[k].gene_id();
					meta_data->gene_names = samples[i]->gene_primary_transcripts[k].gene_name();
					meta_data->protein_ids = samples[i]->gene_primary_transcripts[k].protein_id();
					meta_data->locus_desc = samples[i]->gene_primary_transcripts[k].locus_tag();
					meta_data->description = samples[i]->gene_primary_transcripts[k].description();
					paternalPaternalTest.meta_data = meta_data;
					paternalMaternalTest.meta_data = meta_data;
					maternalPaternalTest.meta_data = meta_data;
					maternalMaternalTest.meta_data = meta_data;
					
#if ENABLE_THREADS
					test_storage_lock.lock();
#endif
					
					pair<SampleDiffs::iterator, bool> paternal_paternal_inserted;
					paternal_paternal_inserted = tests.paternal_paternal_diff_promoter_tests[i][j].insert(make_pair("CurrentPaternal:PreviousPaternal;"+desc,paternalPaternalTest)); 
					paternal_paternal_inserted.first->second = paternalPaternalTest;
					
					pair<SampleDiffs::iterator, bool> paternal_maternal_inserted;
					paternal_maternal_inserted = tests.paternal_maternal_diff_promoter_tests[i][j].insert(make_pair("CurrentPaternal:PreviousMaternal;"+desc,paternalMaternalTest)); 
					paternal_maternal_inserted.first->second = paternalMaternalTest;
					
					pair<SampleDiffs::iterator, bool> maternal_paternal_inserted;
					maternal_paternal_inserted = tests.maternal_paternal_diff_promoter_tests[i][j].insert(make_pair("CurrentMaternal:PreviousPaternal;"+desc,maternalPaternalTest)); 
					maternal_paternal_inserted.first->second = maternalPaternalTest;
					
					pair<SampleDiffs::iterator, bool> maternal_maternal_inserted;
					maternal_maternal_inserted = tests.maternal_maternal_diff_promoter_tests[i][j].insert(make_pair("CurrentMaternal:PreviousMaternal;"+desc,maternalMaternalTest)); 
					maternal_maternal_inserted.first->second = maternalMaternalTest;
					
#if ENABLE_THREADS
					test_storage_lock.unlock();
#endif
				}
				
				// Differential coding sequence output
				for (size_t k = 0; k < samples[i]->gene_cds.size(); ++k)
				{
					const AlleleAbundanceGroup& curr_abundance = filtered_cds_groups[i][k];
					const AlleleAbundanceGroup& prev_abundance = filtered_cds_groups[j][k];
					const string& desc = curr_abundance.description();
					
					bool paternal_paternal_enough_reads = (group_has_record_above_thresh_paternal(curr_abundance) &&
														   group_has_record_above_thresh_paternal(prev_abundance) &&
														   group_has_paternal_record_badly_fit(curr_abundance) == false &&
														   group_has_paternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					bool paternal_maternal_enough_reads = (group_has_record_above_thresh_paternal(curr_abundance) &&
														   group_has_record_above_thresh_maternal(prev_abundance) &&
														   group_has_paternal_record_badly_fit(curr_abundance) == false &&
														   group_has_maternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					bool maternal_paternal_enough_reads = (group_has_record_above_thresh_maternal(curr_abundance) &&
														   group_has_record_above_thresh_paternal(prev_abundance) &&
														   group_has_maternal_record_badly_fit(curr_abundance) == false &&
														   group_has_paternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					bool maternal_maternal_enough_reads = (group_has_record_above_thresh_maternal(curr_abundance) &&
														   group_has_record_above_thresh_maternal(prev_abundance) &&
														   group_has_maternal_record_badly_fit(curr_abundance) == false &&
														   group_has_maternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					
					SampleDifference paternalPaternalTest;
					paternalPaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														paternal_paternal_enough_reads,
														false,false,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference paternalMaternalTest;
					paternalMaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														paternal_maternal_enough_reads,
														true,false,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalPaternalTest;
					maternalPaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														maternal_paternal_enough_reads,
														false,true,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalMaternalTest;
					maternalMaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														maternal_maternal_enough_reads,
														true,true,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					
					// The filtered group might be empty, so let's grab metadata from
					// the unfiltered group
					shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
					
					meta_data->gene_ids = samples[i]->gene_cds[k].gene_id();
					meta_data->gene_names = samples[i]->gene_cds[k].gene_name();
					meta_data->protein_ids = samples[i]->gene_cds[k].protein_id();
					meta_data->locus_desc = samples[i]->gene_cds[k].locus_tag();
					meta_data->description = samples[i]->gene_cds[k].description();
					paternalPaternalTest.meta_data = meta_data;
					paternalMaternalTest.meta_data = meta_data;
					maternalPaternalTest.meta_data = meta_data;
					maternalMaternalTest.meta_data = meta_data;
					
#if ENABLE_THREADS
					test_storage_lock.lock();
#endif
					pair<SampleDiffs::iterator, bool> paternal_paternal_inserted;
					paternal_paternal_inserted = tests.paternal_paternal_diff_cds_tests[i][j].insert(make_pair("CurrentPaternal:PreviousPaternal;"+desc,paternalPaternalTest)); 
					paternal_paternal_inserted.first->second = paternalPaternalTest;
					
					pair<SampleDiffs::iterator, bool> paternal_maternal_inserted;
					paternal_maternal_inserted = tests.paternal_maternal_diff_cds_tests[i][j].insert(make_pair("CurrentPaternal:PreviousMaternal;"+desc,paternalMaternalTest)); 
					paternal_maternal_inserted.first->second = paternalMaternalTest;
					
					pair<SampleDiffs::iterator, bool> maternal_paternal_inserted;
					maternal_paternal_inserted = tests.maternal_paternal_diff_cds_tests[i][j].insert(make_pair("CurrentMaternal:PreviousPaternal;"+desc,maternalPaternalTest)); 
					maternal_paternal_inserted.first->second = maternalPaternalTest;
					
					pair<SampleDiffs::iterator, bool> maternal_maternal_inserted;
					maternal_maternal_inserted = tests.maternal_maternal_diff_cds_tests[i][j].insert(make_pair("CurrentMaternal:PreviousMaternal;"+desc,maternalMaternalTest)); 
					maternal_maternal_inserted.first->second = maternalMaternalTest;
					
#if ENABLE_THREADS
					test_storage_lock.unlock();
#endif
				}
				
				// Differential splicing of primary transcripts
				for (size_t k = 0; k < samples[i]->primary_transcripts.size(); ++k)
				{
					const AlleleAbundanceGroup& curr_abundance = filtered_primary_trans_groups[i][k];
					const AlleleAbundanceGroup& prev_abundance = filtered_primary_trans_groups[j][k];
					const string& desc = curr_abundance.description();
					
					bool paternal_paternal_enough_reads = (group_has_record_above_thresh_paternal(curr_abundance) &&
														   group_has_record_above_thresh_paternal(prev_abundance) &&
														   group_has_paternal_record_badly_fit(curr_abundance) == false &&
														   group_has_paternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					bool paternal_maternal_enough_reads = (group_has_record_above_thresh_paternal(curr_abundance) &&
														   group_has_record_above_thresh_maternal(prev_abundance) &&
														   group_has_paternal_record_badly_fit(curr_abundance) == false &&
														   group_has_maternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					bool maternal_paternal_enough_reads = (group_has_record_above_thresh_maternal(curr_abundance) &&
														   group_has_record_above_thresh_paternal(prev_abundance) &&
														   group_has_maternal_record_badly_fit(curr_abundance) == false &&
														   group_has_paternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					bool maternal_maternal_enough_reads = (group_has_record_above_thresh_maternal(curr_abundance) &&
														   group_has_record_above_thresh_maternal(prev_abundance) &&
														   group_has_maternal_record_badly_fit(curr_abundance) == false &&
														   group_has_maternal_record_badly_fit(prev_abundance) == false &&
														   curr_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
														   prev_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
					
					SampleDifference paternalPaternalTest;
					paternalPaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														paternal_paternal_enough_reads,
														false,false,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference paternalMaternalTest;
					paternalMaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														paternal_maternal_enough_reads,
														true,false,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalPaternalTest;
					maternalPaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														maternal_paternal_enough_reads,
														false,true,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					SampleDifference maternalMaternalTest;
					maternalMaternalTest = get_ds_tests(prev_abundance, 
														curr_abundance,
														maternal_maternal_enough_reads,
														true,true,(curr_abundance.is_allele_informative() && prev_abundance.is_allele_informative()));
					// The filtered group might be empty, so let's grab metadata from
					// the unfiltered group
					shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
					
					meta_data->gene_ids = samples[i]->primary_transcripts[k].gene_id();
					meta_data->gene_names = samples[i]->primary_transcripts[k].gene_name();
					meta_data->protein_ids = samples[i]->primary_transcripts[k].protein_id();
					meta_data->locus_desc = samples[i]->primary_transcripts[k].locus_tag();
					meta_data->description = samples[i]->primary_transcripts[k].description();
					paternalPaternalTest.meta_data = meta_data;
					paternalMaternalTest.meta_data = meta_data;
					maternalPaternalTest.meta_data = meta_data;
					maternalMaternalTest.meta_data = meta_data;
					
#if ENABLE_THREADS
					test_storage_lock.lock();
#endif
					
					pair<SampleDiffs::iterator, bool> paternal_paternal_inserted;
					paternal_paternal_inserted = tests.paternal_paternal_diff_splicing_tests[i][j].insert(make_pair("CurrentPaternal:PreviousPaternal;"+desc,paternalPaternalTest)); 
					paternal_paternal_inserted.first->second = paternalPaternalTest;
					
					pair<SampleDiffs::iterator, bool> paternal_maternal_inserted;
					paternal_maternal_inserted = tests.paternal_maternal_diff_splicing_tests[i][j].insert(make_pair("CurrentPaternal:PreviousMaternal;"+desc,paternalMaternalTest)); 
					paternal_maternal_inserted.first->second = paternalMaternalTest;
					
					pair<SampleDiffs::iterator, bool> maternal_paternal_inserted;
					maternal_paternal_inserted = tests.maternal_paternal_diff_splicing_tests[i][j].insert(make_pair("CurrentMaternal:PreviousPaternal;"+desc,maternalPaternalTest)); 
					maternal_paternal_inserted.first->second = maternalPaternalTest;
					
					pair<SampleDiffs::iterator, bool> maternal_maternal_inserted;
					maternal_maternal_inserted = tests.maternal_maternal_diff_splicing_tests[i][j].insert(make_pair("CurrentMaternal:PreviousMaternal;"+desc,maternalMaternalTest)); 
					maternal_maternal_inserted.first->second = maternalMaternalTest;
					
#if ENABLE_THREADS
					test_storage_lock.unlock();
#endif
				}
			}
		}
	}
	else
	{// do only imprinting tests: pat_vs_mat
		//bool multi_transcript_locus = samples[0]->transcripts.abundances().size() > 1;
			
		//bool enough_reads = (samples[0]->cluster_mass >= min_read_count);
//      
		for (size_t k = 0; k < samples[0]->transcripts.abundances().size(); ++k)
		{
			const Abundance& curr_abundance = *(samples[0]->transcripts.abundances()[k]);
			const string& desc = curr_abundance.description();
			
			FPKMTrackingTable::iterator itr = tracking.isoform_fpkm_tracking.find(desc);
			assert (itr != tracking.isoform_fpkm_tracking.end());
						
			bool curr_paternal_enough_reads = false;
			bool curr_maternal_enough_reads = false;
			
			if ((curr_abundance.num_paternal_fragments()+curr_abundance.num_maternal_fragments()) && (curr_abundance.paternal_effective_length() || curr_abundance.maternal_effective_length()))
			{
				double paternal_frags_per_kb = curr_abundance.num_paternal_fragments() / (curr_abundance.paternal_effective_length() / 1000.0);
				double maternal_frags_per_kb = curr_abundance.num_maternal_fragments() / (curr_abundance.maternal_effective_length() / 1000.0);
				if (paternal_frags_per_kb >= min_read_count)
					curr_paternal_enough_reads = true;
				if (maternal_frags_per_kb >= min_read_count)
					curr_maternal_enough_reads = true;
			}
			if (curr_paternal_enough_reads)
			{
				if (is_badly_fit_paternal(curr_abundance))
					curr_paternal_enough_reads = false;
			}
			if (curr_maternal_enough_reads)
			{
				if (is_badly_fit_maternal(curr_abundance))
					curr_maternal_enough_reads = false;
			}

			SampleDifference paternalMaternalTest;
			paternalMaternalTest = get_de_tests("Paternal:Maternal;"+desc,
												itr->second.paternal_fpkm_series[0], 
												itr->second.maternal_fpkm_series[0],
												(curr_paternal_enough_reads || curr_maternal_enough_reads),
												curr_abundance.is_allele_informative());
			
			
#if ENABLE_THREADS
			test_storage_lock.lock();
#endif
			pair<SampleDiffs::iterator, bool> paternalMaternalInserted; 
			paternalMaternalInserted = tests.parental_isoform_de_tests[0][0].insert(make_pair("Paternal:Maternal;"+desc,
																					 paternalMaternalTest)); 
			
			shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
			
			meta_data->gene_ids = curr_abundance.gene_id();
			meta_data->gene_names = curr_abundance.gene_name();
			meta_data->protein_ids = curr_abundance.protein_id();
			meta_data->locus_desc = curr_abundance.locus_tag();
			meta_data->description = curr_abundance.description();
			paternalMaternalInserted.first->second.meta_data = meta_data;
			
#if ENABLE_THREADS
			test_storage_lock.unlock();
#endif
		}
		
		for (size_t k = 0; k < samples[0]->cds.size(); ++k)
		{
			const Abundance& curr_abundance = samples[0]->cds[k];
			
			const string& desc = curr_abundance.description();
			FPKMTrackingTable::iterator itr = tracking.cds_fpkm_tracking.find(desc);
			assert (itr != tracking.cds_fpkm_tracking.end());
								
			bool curr_paternal_enough_reads = false;
			bool curr_maternal_enough_reads = false;
								
			if ((curr_abundance.num_paternal_fragments()+curr_abundance.num_maternal_fragments()) && (curr_abundance.paternal_effective_length() || curr_abundance.maternal_effective_length()))
			{
				double paternal_frags_per_kb = curr_abundance.num_paternal_fragments() / (curr_abundance.paternal_effective_length() / 1000.0);
				double maternal_frags_per_kb = curr_abundance.num_maternal_fragments() / (curr_abundance.maternal_effective_length() / 1000.0);
				if (paternal_frags_per_kb >= min_read_count)
					curr_paternal_enough_reads = true;
				if (maternal_frags_per_kb >= min_read_count)
					curr_maternal_enough_reads = true;
			}
			
			if (curr_paternal_enough_reads)
			{
				if (is_badly_fit_paternal(curr_abundance))
					curr_paternal_enough_reads = false;
			}
			if (curr_maternal_enough_reads)
			{
				if (is_badly_fit_maternal(curr_abundance))
					curr_maternal_enough_reads = false;
			}
			SampleDifference paternalMaternalTest;
			paternalMaternalTest = get_de_tests("Paternal:Maternal;"+desc,
												itr->second.paternal_fpkm_series[0], 
												itr->second.maternal_fpkm_series[0],
												(curr_paternal_enough_reads || curr_maternal_enough_reads),
												curr_abundance.is_allele_informative());
			
#if ENABLE_THREADS
			test_storage_lock.lock();
#endif
			
			pair<SampleDiffs::iterator, bool> paternalMaternalInserted; 
			paternalMaternalInserted = tests.parental_cds_de_tests[0][0].insert(make_pair("Paternal:Maternal;"+desc,
																				 paternalMaternalTest)); 
			shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
			
			meta_data->gene_ids = curr_abundance.gene_id();
			meta_data->gene_names = curr_abundance.gene_name();
			meta_data->protein_ids = curr_abundance.protein_id();
			meta_data->locus_desc = curr_abundance.locus_tag();
			meta_data->description = curr_abundance.description();
			paternalMaternalInserted.first->second.meta_data = meta_data;
			
#if ENABLE_THREADS
			test_storage_lock.unlock();
#endif
		}
		
		for (size_t k = 0; k < samples[0]->primary_transcripts.size(); ++k)
		{
			const Abundance& curr_abundance = samples[0]->primary_transcripts[k];
			
			const string& desc = curr_abundance.description();
			FPKMTrackingTable::iterator itr = tracking.tss_group_fpkm_tracking.find(desc);
			assert (itr != tracking.tss_group_fpkm_tracking.end());
						
			bool curr_paternal_enough_reads = false;
			bool curr_maternal_enough_reads = false;
			
			if ((curr_abundance.num_paternal_fragments()+curr_abundance.num_maternal_fragments()) && (curr_abundance.paternal_effective_length() || curr_abundance.maternal_effective_length()))
			{
				double paternal_frags_per_kb = curr_abundance.num_paternal_fragments() / (curr_abundance.paternal_effective_length() / 1000.0);
				double maternal_frags_per_kb = curr_abundance.num_maternal_fragments() / (curr_abundance.maternal_effective_length() / 1000.0);
				if (paternal_frags_per_kb >= min_read_count)
					curr_paternal_enough_reads = true;
				if (maternal_frags_per_kb >= min_read_count)
					curr_maternal_enough_reads = true;
			}
			
			if (curr_paternal_enough_reads)
			{
				if (is_badly_fit_paternal(curr_abundance))
					curr_paternal_enough_reads = false;
			}
			if (curr_maternal_enough_reads)
			{
				if (is_badly_fit_maternal(curr_abundance))
					curr_maternal_enough_reads = false;
			}
					
			SampleDifference paternalMaternalTest;
			paternalMaternalTest = get_de_tests("Paternal:Maternal;"+desc,
												itr->second.paternal_fpkm_series[0], 
												itr->second.maternal_fpkm_series[0],
												//tests.tss_de_tests[0][0],
												(curr_paternal_enough_reads || curr_maternal_enough_reads),
												curr_abundance.is_allele_informative());
			
			
#if ENABLE_THREADS
			test_storage_lock.lock();
#endif
			pair<SampleDiffs::iterator, bool> paternalMaternalInserted; 
			paternalMaternalInserted = tests.parental_tss_group_de_tests[0][0].insert(make_pair("Paternal:Maternal;"+desc,
																				 paternalMaternalTest)); 
			
			shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
			
			meta_data->gene_ids = curr_abundance.gene_id();
			meta_data->gene_names = curr_abundance.gene_name();
			meta_data->protein_ids = curr_abundance.protein_id();
			meta_data->locus_desc = curr_abundance.locus_tag();
			meta_data->description = curr_abundance.description();
			paternalMaternalInserted.first->second.meta_data = meta_data;
			
#if ENABLE_THREADS
			test_storage_lock.unlock();
#endif
		}

		for (size_t k = 0; k < samples[0]->genes.size(); ++k)
		{
			const AlleleAbundanceGroup& curr_abundance = samples[0]->genes[k];
			const string& desc = curr_abundance.description();
			FPKMTrackingTable::iterator itr = tracking.gene_fpkm_tracking.find(desc);
			assert (itr != tracking.gene_fpkm_tracking.end());
						
			bool curr_paternal_enough_reads = false;
			bool curr_maternal_enough_reads = false;
			
			if ((curr_abundance.num_paternal_fragments()+curr_abundance.num_maternal_fragments()) && (curr_abundance.paternal_effective_length() || curr_abundance.maternal_effective_length()))
			{
				double paternal_frags_per_kb = curr_abundance.num_paternal_fragments() / (curr_abundance.paternal_effective_length() / 1000.0);
				double maternal_frags_per_kb = curr_abundance.num_maternal_fragments() / (curr_abundance.maternal_effective_length() / 1000.0);
				if (paternal_frags_per_kb >= min_read_count)
					curr_paternal_enough_reads = true;
				if (maternal_frags_per_kb >= min_read_count)
					curr_maternal_enough_reads = true;
			}
			
			SampleDifference paternalMaternalTest;
			paternalMaternalTest = get_de_tests("Paternal:Maternal;"+desc,
												itr->second.paternal_fpkm_series[0], 
												itr->second.maternal_fpkm_series[0],
												//tests.gene_de_tests[0][0],
												(curr_paternal_enough_reads || curr_maternal_enough_reads),
												curr_abundance.is_allele_informative());

#if ENABLE_THREADS
			test_storage_lock.lock();
#endif
			pair<SampleDiffs::iterator, bool> paternalMaternalInserted; 
			paternalMaternalInserted = tests.parental_gene_de_tests[0][0].insert(make_pair("Paternal:Maternal;"+desc,
																				  paternalMaternalTest)); 
			
			shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
			
			meta_data->gene_ids = curr_abundance.gene_id();
			meta_data->gene_names = curr_abundance.gene_name();
			meta_data->protein_ids = curr_abundance.protein_id();
			meta_data->locus_desc = curr_abundance.locus_tag();
			meta_data->description = curr_abundance.description();
			paternalMaternalInserted.first->second.meta_data = meta_data;
			
#if ENABLE_THREADS
			test_storage_lock.unlock();
#endif
		}

		// Skip all the JS based testing for genes with an isoform switch?
		if (!no_js_tests)
		{

			// FIXME: the code below will not properly test for differential
			// splicing/promoter use when a gene (e.g.) occupies two
			// disjoint bundles.  We need to store the covariance matrices (etc)
			// in the FPKMContexts to handle that case properly.
			
			// Differential promoter use
			for (size_t k = 0; k < samples[0]->gene_primary_transcripts.size(); ++k)
			{
				const AlleleAbundanceGroup& curr_abundance = filtered_promoter_groups[0][k];
				const string& desc = curr_abundance.description();
				
				bool paternal_maternal_enough_reads = (group_has_record_above_thresh_paternal(curr_abundance) &&
													   group_has_record_above_thresh_maternal(curr_abundance) &&
													   group_has_paternal_record_badly_fit(curr_abundance) == false &&
													   group_has_maternal_record_badly_fit(curr_abundance) == false &&
													   curr_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
													   curr_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
				
				SampleDifference paternalMaternalTest;
				paternalMaternalTest = get_ds_tests(curr_abundance, 
													curr_abundance,
													paternal_maternal_enough_reads,
													false,true,curr_abundance.is_allele_informative());
				
				// The filtered group might be empty, so let's grab metadata from
				// the unfiltered group
				shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
				
				meta_data->gene_ids = samples[0]->gene_primary_transcripts[k].gene_id();
				meta_data->gene_names = samples[0]->gene_primary_transcripts[k].gene_name();
				meta_data->protein_ids = samples[0]->gene_primary_transcripts[k].protein_id();
				meta_data->locus_desc = samples[0]->gene_primary_transcripts[k].locus_tag();
				meta_data->description = samples[0]->gene_primary_transcripts[k].description();
				paternalMaternalTest.meta_data = meta_data;
				
				
#if ENABLE_THREADS
				test_storage_lock.lock();
#endif
				
				
				pair<SampleDiffs::iterator, bool> paternal_maternal_inserted;
				paternal_maternal_inserted = tests.parental_diff_promoter_tests[0][0].insert(make_pair("Paternal:Maternal;"+desc,paternalMaternalTest)); 
				paternal_maternal_inserted.first->second = paternalMaternalTest;
				
				
#if ENABLE_THREADS
				test_storage_lock.unlock();
#endif
			}
			
			// Differential coding sequence output
			for (size_t k = 0; k < samples[0]->gene_cds.size(); ++k)
			{
				const AlleleAbundanceGroup& curr_abundance = filtered_cds_groups[0][k];
				const string& desc = curr_abundance.description();
				
				
				bool paternal_maternal_enough_reads = (group_has_record_above_thresh_paternal(curr_abundance) &&
													   group_has_record_above_thresh_maternal(curr_abundance) &&
													   group_has_paternal_record_badly_fit(curr_abundance) == false &&
													   group_has_maternal_record_badly_fit(curr_abundance) == false &&
													   curr_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
													   curr_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
				
				
				SampleDifference paternalMaternalTest;
				paternalMaternalTest = get_ds_tests(curr_abundance, 
													curr_abundance,
													paternal_maternal_enough_reads,
													false,true,curr_abundance.is_allele_informative());
				
				// The filtered group might be empty, so let's grab metadata from
				// the unfiltered group
				shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
				
				meta_data->gene_ids = samples[0]->gene_cds[k].gene_id();
				meta_data->gene_names = samples[0]->gene_cds[k].gene_name();
				meta_data->protein_ids = samples[0]->gene_cds[k].protein_id();
				meta_data->locus_desc = samples[0]->gene_cds[k].locus_tag();
				meta_data->description = samples[0]->gene_cds[k].description();
				paternalMaternalTest.meta_data = meta_data;
				
#if ENABLE_THREADS
				test_storage_lock.lock();
#endif
				
				pair<SampleDiffs::iterator, bool> paternal_maternal_inserted;
				paternal_maternal_inserted = tests.parental_diff_cds_tests[0][0].insert(make_pair("Paternal:Maternal;"+desc,paternalMaternalTest)); 
				paternal_maternal_inserted.first->second = paternalMaternalTest;
				
#if ENABLE_THREADS
				test_storage_lock.unlock();
#endif
			}
			
			// Differential splicing of primary transcripts
			for (size_t k = 0; k < samples[0]->primary_transcripts.size(); ++k)
			{
				const AlleleAbundanceGroup& curr_abundance = filtered_primary_trans_groups[0][k];
				const string& desc = curr_abundance.description();
				
				
				bool paternal_maternal_enough_reads = (group_has_record_above_thresh_paternal(curr_abundance) &&
													   group_has_record_above_thresh_maternal(curr_abundance) &&
													   group_has_paternal_record_badly_fit(curr_abundance) == false &&
													   group_has_maternal_record_badly_fit(curr_abundance) == false &&
													   curr_abundance.paternal_FPKM_by_replicate().size() >= min_reps_for_js_test && 
													   curr_abundance.maternal_FPKM_by_replicate().size() >= min_reps_for_js_test);
				
				
				SampleDifference paternalMaternalTest;
				paternalMaternalTest = get_ds_tests(curr_abundance, 
													curr_abundance,
													paternal_maternal_enough_reads,
													false,true,curr_abundance.is_allele_informative());
				
				// The filtered group might be empty, so let's grab metadata from
				// the unfiltered group
				shared_ptr<SampleDifferenceMetaData> meta_data = get_metadata(desc);
				
				meta_data->gene_ids = samples[0]->primary_transcripts[k].gene_id();
				meta_data->gene_names = samples[0]->primary_transcripts[k].gene_name();
				meta_data->protein_ids = samples[0]->primary_transcripts[k].protein_id();
				meta_data->locus_desc = samples[0]->primary_transcripts[k].locus_tag();
				meta_data->description = samples[0]->primary_transcripts[k].description();
				paternalMaternalTest.meta_data = meta_data;
				
#if ENABLE_THREADS
				test_storage_lock.lock();
#endif
				
				pair<SampleDiffs::iterator, bool> paternal_maternal_inserted;
				paternal_maternal_inserted = tests.parental_diff_splicing_tests[0][0].insert(make_pair("CurrentPaternal:CurrentMaternal;"+desc,paternalMaternalTest)); 
				paternal_maternal_inserted.first->second = paternalMaternalTest;
				
				
#if ENABLE_THREADS
				test_storage_lock.unlock();
#endif
			}
		}
	}
}
