#import "@preview/cetz:0.4.2": canvas, draw
#import draw: arc, circle, content, line, rect
#let arrow-style = (
  mark: (end: "stealth", fill: black, scale: 0.7),
  stroke: 0.8pt,
)

=== TODO
 Sparse firing model
 Check the other method (Matlab).
 
  What smoothing to use?

  Process older cells (cheng ren)

 Issue with unstable SIC values.

== 2026g0106

Still struggling to get something reasonable out of the joint regression. I get some reasonble numbers for the proportion of cells with locationg tuning ($~35%$), but very few cells with gaze tuning. It seems like the main trouble with gaze tuning is that the number of effective degrees of freedom is not really affect much by the smoothing. This is in turn is caused by the fact that the offset term is so much larger than the other coefficients in the fit. The result of that is that the term $X W X'$ completely dominates the smoothing term $alpha L$, and so the smoothing does not affect the number of degrees of freedom much at all. One possibility is to explore larger $a$ values. Currently, I'm only going up to an $alpha$ value of 0.01, but some cells do appear to show a continuing improvement in log-likelihood even up to this values. For these cells, it is possible that increasing the smoothing factor will in turn decrease the value of the constant factor, and so reduce the number of degrees of freedom to something a bit more reasonable.

Technically, we should check for convergence of the model, but presumably if it has not converged, it should perform better than the null model. 

== 20260107

I found another bug, or at least an issue, where I was not using the same training and testing indices for each of the variable combinations. This meant that I could not do a paired test between the log-likelihoods for different variable fits, which is necessary if I want to make any statement about whether a cell is more location or more gaze sensitive, for instance. I'm doing another run now where I'm using the original training indices from the gaze objects to re-fit the place, head direction and all the other combinations.

I am also again trying to get a sense of how many cells are modulated by place. I'm going to do this in two mways. The first is the way that the Mao Dun paper did it, namely to compare the log-likelihoods of the test data under the fitted model This avoids (at least to some extent) the issue of overfitting that occurs when introducing more variables (i.e. when looking at combinations of variables), because overfitting on the training data will lead to lower likelihood of the testing data. 
The other approach is to use the likelihood ratio test on the fits to the training data, using an esitmate of the effective number of degrees of freedom to determine significance.

Again, if I just use the log-likelihood of the testing data as a measure, most of the cells are sensitive to place. In fact, for 230 of the cells, across all 10 cross-validation runs, the log-likelihood of the mode was already larger than that of the null model. 

Another bug. The default refinement for gaze was set to 3, rather than 2 which is the one that I'm currently using (since it makes the number of gaze bins comparable to the number of place bins). This meant that when I was fitting models to place responses, the JointOccupancy object I was using had more gaze bins than the one used when fitting to gaze responses. The result of this was that I had more bins to fit for place than for gaze, and so it was not possible to match the training and testing sets. I have no fixed this such that the default refinement for gaze is 2. That means that I need to redo the fits for place, yet again. I also want to try and increase the number of runs from 10 to 20, since for a few cells where the log-likelihood for place/gaz minutely larger than the null-likelihood, I saw more reasonable comparison values when using 20 runs. I might not need to redo the full cross-validation, though. I could just do a single set of runs for the best alpha from the 10 run set.

=== TODO
- #sym.ballot This is pending
- #sym.ballot.check This is done

I'm seeing a persistent issue with number of bins being different from location and gaze, even after changing the default refinement for gaze when loading JointOccupancy objects. Found the bug, I had forgotten to change the defaults in `process_refinements`.

== 20260108

Still waiting for all the cells to complete, but at least I now have 4 computers (work22, work23, work24, work30) running their own group of cells. I could also add work28, of course.

== 20260109

I'm still confounded by the results I'm getting. Basically, way too many cells show up as significantly coding for either place or gaze. A couple of things I want to try out

- #sym.ballot Instead of using a paired test, just compare the overall distribution of likelihood for the null model and the full model
- #sym.ballot.check Check what happens if I shuffle the relationship between spiking activity and behavioural variables. In theory, the likelihood for both models should be identical here.
  - It looks like this is still larger than the null, suggesting that the model is specializising on the training data


Looking at some of the plots, I notice that the log-likelihood of the training data and the testing data tend to be anti-correlated. This means that, whenever the model fits the training data well, it fails to explain the testing data, and vice versa. This again suggests that the model is not generalizing, which could explain the results that I'm seeing. I'm currently checking how prevalent this is. I know that I've seen some cells for which the log-likelihood for training and testing data were positively correlated.

In the meantime, what could cause the specialization? One possibility is lack of mixing between training and testing. Could it also be that insufficient smoothing naturally leads to more specialization? In principle, yes. Right now I'm determining the optimal smoothing factor by find the largest log-likelihood on the test set. If a smaller smoothing factor leads to higher specialization, then that should also be reflected in the log-likelihood; more specialization on the training set should lead to a worse fit on the testing set.  If both the null model and the full model exhibit negative correlations between training and testing goodness-of-fit, it becomes a contest of which one is less negatively correlated. Here, since the full model has more parameters, it should be the most negatively correlated, and hence, it should again exhibit worse goodness-of-fit than the null model. Why does this not happen?

Could it be that the temporal scale is too fine? There are a lot more bins with zero spikes than bins with 1 spike or more. Could this naturally lead to a model that only tries to fit the zeros? To fix this, I can try re-sampling. Currently, I am using the actual sampling of the eyelink data, i.e. 1 ms sampling, and only grouping bins that are identical (i.e. if two identical bins are visited twice, they get grouped together with a duration of 2 ms). For coarses binning, I can do it "standard way", i.e. for each trial, using a fixed window size and compute the average behavioural response (location or gaze) within that window.

One caveat here; even though I am sampling very finely, the mean bin size is around 40ms, i.e. not very different from the 20ms normally used. 

For example, I have one cell for the which total number of spikes is 657 over a period of 1350 s, which amounts to a mean firing rate of 0.5 Hz. This is not an unreasonable number, I think.  Wait, did I just over-compensate for the time window? When I compute $lambda = exp(beta ^ X)$, that #emph[is] the mean spike count in that bin. It doesn't matter what the time window is. In my 

== 20260110

I found a bug in the code where I scaled the estimated firing rate $lambda = exp(beta X)$ by the window size $d t$. This is what caused the log-likelihood for the full model to always be larger than for the null model. Now, it seems, nothing is significant anymore (which I think is slightly better than everything being significant). If the code is bug free (a huge IF), it could be that the imbalance between the number of spike counts that are zero vs non-zero is throwing the model fit off. In other words, it is only fitting the zero counts while ignoring the non-zero counts. One way to reduce the number of zero counts is to be stricter in what bins to include. For instance, by increasing the speed threshold from 1 to 2. This would exclude bins in which the animal was navigating very slowly, in which there also might not be that many spikes. 
(I should also check what bins are actually included, just to verify that I'm only including bins that are visited. )

Another observation; the offset is currently not penalized (I'm setting the corresponding entry in the L matrix to 0), which I think causes the model to favour the offset in explaining the data. I'm goign to try to add a 1 instead, which means I'm using ridge regression. It is possible that the penality should be separate for this term, but for now I'm just going to use the same $alpha$ value as for the rest of the parameters.

== 20260111

This is the most promising cell so far (this is `place_and_view_selective_cells[7]`). 

#figure(image("../figures/hardcastle_place_tuning_p20181102s01a02g045c02.png"), caption: [A cell which was previously identified as having both place and view modulated responses. ])

As we can see in C), the log-likelihood for the full model (y-axis) is only slightly smaller than the log-likelihood for the null model. In addition, it looks like the improvement in log-likelihood does not peak at a smoothing factor of 0.01 (B), but rather continues. Thus, I'm running this cell again, trying smoothing factors up to 1. 
We do see a fairly distinct place-like field on the inner side of the top-left pillar, which is also promising.

== 20260112

I also see some promising results for cell `place_selective_cells[end-8]`, show below

#figure(image("../figures/hardcastle_place_tuning_p20181026s01a02g035c04.png"), caption: [A cell with a previously identified place field.])

This cell is still not significant, though, so I'm checking whether increasing the speed threshold changes anything. 

I still think that the number of bins with zero spikes might be an issue. For instance, for the above cell, there are a total of 31285 bins with no spikes, 26 bins with 1 spike and 2 bins with 2 spikes. This is after restricting the bins we use to only those actually visited, and with sufficiently higher (i.e. $>1$) speed. This does mean that this cell has a really low mean firing rate of only 0.02 Hz.

Perhaps we simply can't model such sparse cells with a Poisson distribution. Perhaps try the zero inflated Poisson model for these cells?

I recognize that perhaps it is not fruitful to spend too much time on this. My issue now is that the Hardcastle method doesn't seem ...

Now the problem seems to be that the offset is weighted much less than the other parameters.

One possible issue is that I've been using the unnormalized laplacian. This should suppress the coefficients, though, rather than enhance them, since it puts a higher cost on the invidual parameters compared to the intercept.  Channging to a normalized laplacian doesn't appear to make much difference.

Is it possible that the spatial mapping is wrong? Currently running a basic experiment where I manually set the spike counts for some subset of the bins and the run the GLM fit on that. This should produce a map replicating that pattern. If it doesn't, then something is clearly wrong the with GLM procedure itself. If it does, then then next stop is to verity that the spatial mapping actually works, i.e. that the $X$ matrix that goes into the fit procedure actually represents the true spatial variables (I have of course checked this at various stages already and everything seems to check out, but I guess I should create a rigorous test set to actually prove it.) 

Phew, it does work.

== 20260113

#figure(image("../figures/path_placebin_correspondence_p20181102s01_trial2.png"), caption: [Illustration of correspoondence between spatial binning and actual path. The black dots mark the neareat spatial bin (grid) to the actual path (red).])

As the above figure shows, the mapping from actual spatial path (red) to spatial binning (black dots) also makes sense. This was just a sanity check, basically, to make sure things do not get messed up in the binning.

So, to conclude, we know that the Hardcastle method can uncover the spatial firing pattern, and we know that the spatial binning preserves the original spatial relationship between spikes and location. Why, then, does the method not work for actual cells?

The general problem I seem to have it that log-likelihood of the null model is always higher than for the full model This suggests that the null model actually generalizes better than the full model, which makes sense if the extra degrees of freedom in the full model are mainly just fitting noise in the training data. In general, applying smoothing should alleviate this. Thus, with more smoothing, the model should generalize better. 

But weirdly, it appears that the in-sample vs out-of-sample log-likelihood become more anti-correlated with increasing smoothing. 

Something weird might be going on with the normalization as well; if I normalize the laplacian, the results lookk very strange.

Actually normalization removes the distinction between edge and internal vertices, it seems, so think I'm just not going to normalize. 

For comparsion to the null model, perhaps I can try applying the same smoothing. Actually the penalty has very little effect on the estimated firing rate. For reference, with a penalty term on the square of the firing rate $lambda$, the MLE estimate is

$ lambda = frac(-T + sqrt(T^2 + 8 alpha N),4 alpha)) $

where 

$T = sum_i Delta t_i $ and $N=sum_i n_i$ 

== 20260114

I decided to use the same fitting framework for the null model. Basicaly, what I did above is incorrect because the penalty is not on the firing rate itself, but rather on the $beta$ coefficient, which is the log of the firing rate. I'm not sure if it is possible to find a closed-form expression for $lambda$ under such constraints, so I decided to just train a model exactly like I do for the location and gaze variables, but using only offset. I think this should be the most fair comparison. Indeed, for the cell I was looking at earlier, it did end up being significantly modulated by space, i.e. the log-likelihood of the full model did exceed that of the null model for the testing set. It remains to be seen whether this will now again always be true, i..e that all cells are now significant, in which case I am back to square one. I'm going to run a couple of more cells from the place selective set with decent firing rates, and then try a couple of cells that were not deemed place selective from our previous analysis. Basically, I want to see some diversity, where some cells are significant and some are not. Once that works, but plane is to, yet again run this analysis for all the cells.  I also decided to decrease the resolution a bit to speed up the computation. That means that I'm using a refinement of 2 instead of 3 for place, which means that I'm using 336 bins for the floor instead of 1344 bins. After the preliminary run, I can redo this at higher resolution if needed.

Actually, we should be taking the time window into account as well, since I am not using the same counting window for each bin. To do this, I change the defining of $lambda$ to $lambda = delta t_i exp(beta X_(:,i))$ where $delta t_i$ is the time window for bin $i$. 


== 20260115

I still haven't a cell that does not show significant spatial modulation. This is troubling. Welll, it could be because I am using the smoothing weight from the full model in the null model. I think it makes more sense to find a different smoothing parameter for the null model via the same cross-validation as for the full model. That way, I'm not unduly constraining the null model; after all, they are different models; in the null model I'm trying to find a (regularized) constant firing rate that best explains the spiking activity, while in the full model, I'm trying to find the best $beta$ coefficients that map the behavioural variables (i.e. spatial location or gaze) to a per-bin firing rate. 

=== Spatial representation

I think it would make sense to classify a cell as spatially selective if the SIC for the raw firing rates exceeds that of the shuffle. This is apparently what Daigle et al (2025) did. I think it makes sense, as the dependency on space should not depend on the smoothing; if a cell is spatially tuned, the spikes should reflect that in the form of an excess of spatial orgnaization compared to a shuffled representation. Smoothing should then be applied to those cells that turn out be spatially tuned to quantity that cell's spatial field. 


== 20260117

In their preprtint, Daigle et al (2025) show that most cells with place fields exhibit them in task relevant locations, i..e they are not randomly distrbuted in space. I think we see something similar, but it would be good to check whether place fields are concentrated around the poster locations in our case.

== 20260119

It looks like around 40-50% of the cells are identified as place cells using the SIC computatation for the raw firing rates. This is much higher, obiviously, than the roughly 10% that we identified before using SIC on the smoothed firing rates.

== 20260120

The new plan is to try and run the simpler spatial tuning method first. That is, compute SIC for the raw maps, establish significance via shuffling, and then label those cells for which the SIC is significant as spatially tuned. We can do that for gaze tuning as well, and also for head direction tuning. Then, for those cells that are deemed spatially tuned, identify their response fields by find bins where the firing rate is significantly higher than for for the shuffled surrogates. This will most likely results in many discontiguous fields, so we need to find a awy to consolidate them. 

== 20260121

There are at least two ways to approach the place responses. We can look at significant bins as we coarsen the binning of the arena; we can then look at "lfetimes" of bin signifiance, akin to what one does in conformal homology. There could be bins that are significant at a very high resolution, but that lose their selectivity once we group that bin's responses with those of neighbouring bins. Bins with long lifetimes would be those that remain part of a significant cluster of bins as the scale coarsens. 
Another way is to increase the level of smoothing and then find the bins that survive across the longest range of smoothing parameters. 
Do we also require the SIC to be significant? In a way the SIC avoid the issue of multiple comparisons, since if we are only looking at single bins being significant, it matters that we are actually performing one statistical test per bin, and so we need to do corrections for these multiple comparisons. If we just do SIC, then we are essentially asking whether there is spatial information in the spiking above that of the shuffled surrogates. We can have isolated bins that are significant, but the spatial information they represent might not be significant. 

How do we define the "true" place field though? If we have individual bins with firing rates higher than for the shuffled, but the overall distribution of firing rates does not convey significant information about space, I think we do not way to label that cell as a place cell, at least not at the scale. However, if for a different scale, the firing rates do convery information about space, we can then examine the significant bins. 

Does it make sense to choose the largest scale at which we do see significant spatial information and define place fields from that? What about place fields that only exist at a smaller scale? For instance, if at a resolution of 1 we find that the overall map contains spatial information, and that there a number of significant bins, but that one (or a group of) those bins disappears if we decrease the scale (with the assumption that we still have spatial information at that scale)? I think it would be legit to say that those bins constitute smaller place fields.

== 20260123

Is it perhaps the case that using the coarsest scale we always get a significant SIC score? Mainly because there are only so many ways in which to arrange the bins. Well, we could still have a situation where shuffling basically just results in all the bins being equally occupied. 

== 20260125

Trade-off; using a coarses scale seems to substantially increase the number of cells deemed to be spatially selective. This larger scale of course means that we could lose precision; we won't detect smaller fields. By looking across scales, this problem is largely fixed. Perhaps, then, rather than using one scale to report the number of cells with spatial selectivity, we should instead just reporrt the number of place fields; units without spatial selectivity would simply have zero fields. The problem with this is that, since we have to look across multiple spatial scales for each cell, this analysis takes a lot of time.

== 20260129

I might change to laplace smoothing, mainly because it is more efficient than gaussian smoothing; it only requires matrix multiplications, rather than finding all pairwise distances between elements. With proper tweaking of parameters, I think I can get it to replicate gaussian smoothing pretty closely. At the end of the day, the precise method of smoothing doesn't matter; we only want a somewhat continuous representation of a a cell's firing rate as a function of spatial location and gaze position.

Regarding computing SIC on the joint place and view responses, that appears to be taking forever.

== 20260130

For my test cell, `allcelldirs[39]`, increasing the amount of laplace smoothing, going from 1000 iteratations to 2000 iterations with $alpha$ fixed at 0.1, results in the SIC for the spatial map going from being significant to being non-significant

#figure(image("../figures/p20181031s01a01g030c01_place_map_increase_laplace_smoothing.png"),
       caption: [Raw place field (left) and place field with increasing degree of Laplace smoothing. The bottom rows shows the distribution of SIC scores for the map shown on top. When increasing the nubmer of iterations to 2000, this cell no longer shows a significant difference between the data (red dot) and the surrogates.])


== 20260205

What do we do with cells which are spatially selective according to the SIC, but for which no field is identified? 

How do we think about the conjunctions? If the fields that we observe are really the reults of a conjunction, i.e. if firing rate is modulated by a specific combination of place and view input, then if we condition on the place field, we should see higher SIC in the view responses than if we do not condition? We need some kind of statistics for that, though. If we look at the difference in SIC, we can do a shuffle test, where we basically randomly permute the label (infield vs outfield) and recompute the SIC. That is, we are mixing responses from in-field and out-field. Actually what makes more sense is to grab a random number of bins outside the place field and compare the SIC for gaze when conditioning on those bins

== 20260209

I'm a bit puzzled that we only got two cell with sinificant conjunctions, and those conjunctions looked really strange. Now I'm wondering whether using SIC determine the conjunction is not really the right approach. What we are really asking is whether the observed place and view activity represent a joint response, or if they are simply separate place and view activity maps. If they were conjunctions, then if we conditioned on the place may, say, the activity in the view field should be strenghtened compared to when we condition on other place bins. That is because, in the extreme case, the view field is only active when the the subject is standing in particular location, and so if we condition on other locations, the view field should be weakened, or in even absent. Thus, a more straight-forward test is to simply compute the distribution of firing rates in the identified view field when conditioning on a place field, or when conditioning on other bins.


== 20260211

My priority right now is to try and figure out why we got so few place selective cells in our previous analysis. There, we were using adaptive smoothing quite aggressively. I'm not able to reproduce the result using the same approach, though there are a lot of differences in the approaches. Perhaps a more fruitful way forward is to just make sure that what I'm currently getting is real. 

What could inflate the results? 
One possibility is the random shifts. I'm currently randomly shifting the entire spike train.

=== Conjunction
So far, I'm not seeing a lot of conjunctions. That is, if I comare the mean firing rate within a view field when conditioning on a place field, compared to when conditioning on a location not in the place field, I find that only 5 (of 71 with selectivity for both place and view) cells show a significant difference for at least one view field. This, again, is in contrast to what we found before, we found that most of the cells with both place and view activity were conjunctive. 
The method we used earlier was a bit different, though. I think we just compared the firing rates when conditioning on place, but not specifically in the view field. That is, is the median firing rate in the view space higher when we condition on a place field versus if we don't.  This has the advantage of not requiring a defined place field; it is simply asking whether the overall distribution of firing rates in the view space is influenced by the firing rate in location space. Note that it does require there to at least one place field so that we have something to condition on. 

We could try an ANOVA type analysis for the conjunction as well. Though it wouldn't be exactly like ANOVA. One possibility is to treat the firing rate maps (place and view, and perhaps even head direction) as the rate in a Poisson process. That is, we count the number of spikes and check whether a Poisson model with firing rates describe the firing better. But what if we just used ANOVA? That is, for each spike count entry, we have a place bin and a view bin, as well as a head direction bin. We can then do a linear model to see which variable influences the firing rate the most.

Mistral just suggested computing mutual information between place and view spaces, which might be interesting, though this seems similar to the SIC analysis we had in mind earlier. That SIC analysis took a while to compute. What we could do, though, is to use nearest neighbour approximation to compute the entropy.

== 20260212

 I have made some progress in my effort to try and classify the dependency of place and view in determining the firing rate. Basically, I want to try and use the mutual information, computed according to the KSG estimator, as a meausre of dependence. The advantage of this method is that it allows a directy comparison of the distribution of firing rates when conditioning on both place and view, and the firing rates conditioned only on place, or only on view. My current issue is that I would like to create a simple model to validate that the estimator actually works. I need a way to generate 1D variables with varying dependency on some underyling variables. For instance, the variable could only depend on X, only on Y, or on some combination of the two. I guess I can just use a linear model; generate gaussian variables X and Y and compute a third variable Z as a combination of X and Y.

Hm, I think it doesn't make sense to think of the dependency like this. We are not comparing rates; we are comparing rates conditioned on different things. Perhaps it comes back to computing SIC; Is the SIC bigger when conditioning on different things.


== 20260219

We have essentially two levels of analysis for place/view fields. The first is whether a cell is overall selective to space/view. For this we use SIC on smoothed maps and compared the true SIC value to that of shifted surrogates. The next level is to find the actual place/view fields. What we did in the past was basically to use the same smoothed maps and estimate the outline of significant peaks. This involves choosing value for what amplitudes should be considered peaks, and how large an area to grow around each peak.
An alternative method is to compare the firing rates of each bin the smoothed map with firing rates in that bin for the surrogates. Bins for which the true firing rate exceeds some threshold determined by the firing rate of the surrogates are considered significant. We then group adjacent bins and (tentatively) call these the place/view fields. We could also do a secondary, cluster type analysis, to see who many adjacent bins should be considered a significant grouping. A single bin does probably not consitute a field, but what about two adjacent bins? The approach here would be to randomly draw $n$ siginficant bins from a total of $N$ possible bins and see how many times we get clusters of size $n_s$.

== 20260220

Can we have place selectivity without place fields? Meaning, if a cell is marked as selective to space via SIC, but then has no bins where the firing rate xceeds the 99th percentile of the surrogates, do we simply call this a false positive? Actually, a cell can be spatially selective without being a place cell, I think. 


== 20260224

Why do the field locations not always coincide with the peak of the smoothed map? Actually, they do. There was a bug in the plotting code, where I was using gaussian smoothing for the display, but laplace smoothing for the actual field calculations.


== 20260226

Going to back to the UnityRaytraceData object to figure out how we can segment it into sessions. It looks like the timestamps are reset between every trial

== 20260227

It seems we do get a lot fewer conjunction cells in the new analysis than what we got before. In the previous run, most of the cells with both place and view responses were also conjunctive, whereas currently, only about 5.5% of the cells (6/109) have significant place-conditioned view responses and 11.1% (12/109) have significant view-conditioned place responses. Are we just being more stricter now? We do get a lot more both place and view cells, and therefore also more cells with both place and view responses. This is how I currently compute the conjunctions.
If we are computing view-conditioned place responses, I first find all the view fields of the cell. A view field is defined as a group of bins where the firing rate in the smoothed map exceeds the 99th percentile of the corresponding firing rate in the surrogates for that same bin. In addition, to avoid spurious small fields, I also require that a putative grouping must be larger than what we get by chance by randomly distributing the significant bins among the total number of bins. 
Once we identified the view fields, we now compute the place map while conditioning on those fields. That is, we compute the place map using those place bins for which the view bins fell into one of the view fields. We do this separately for each view bin. As a comparison, we compute the place map if we conditioned on view bins that do not fall in any of the identified view fields. We then look at the average firing rate within each of the original place fields, for the view conditioned responses, and compare it to the distribution of average firing rates in the same place field when conditioning on view bins outisde the view fields. If the view-conditioned firing rate was larger than the 99th percentile of the randomly conditioned firing rates, we identity the cell as having a significant view conditioned place response. We followed an analogous procedure for identifying significant place-conditioned view responses. 
This procedure specifically identifies whether a place (or view) field is the result of a conjunction between view (place) and place (view). However it is also possible that a cell's view conditioned place response is in general different than the randomly conditioned response. Thus, we also ran an analysis where we computed the SIC for the view conditioned response and compared it to the randomly conditioned response. 

Do we do SIC on the raw or the smoothed responses? For place and view fields, we always use the smoothed responses. We could do the joint smoothing of place and view. This essentially gives us a smooth estimate of the firing rate given both place and view. We can then condition this on particular parts of space or view.

$ lambda(p|v) = sum_v lambda(p,v) p(v) $

What is $p(v)$, i.e. the probability of a given view bin $v$? It could just be relative amount time spent in that bin, i.e.

$ p(v) = frac(sum_i (w_i=v),sum(w_i)) $

where $w_i$ is the time spent in bin $i$.

== 20260302

What are some potential issues with the approach to identifying conjunctions? 
- Field identification
  - This is obviously critical, since the fields determine how we condition. If the fields are not well defined, the conditioning is not meaningful, and we might not be able to identify difference between the in-field and out-field responses. 
  - The shuffling is somehow not being done correctly
    - If the shuffling ends up mixing responses associated with e.g. place field with those that are not, then we do not expect to see differences.
- Using raw instead of smoothed responses. 
  - The raw responses could be subject to noise, and we could have noisy spikes both in the place-conditioned and the unconditioned responses. The smoothed responses could be cleaner. The downside of this is that it takes relatively long to smooth the joint responses.


  == 20260304

  Something strange about speed processing. When doing space along, with has lower temporal resolution, we can just filter by speed directly; points that are stationary, e.g. when the subject is not moving, can be excluded. For combined view and place, though, this doesn't really make sense. We do want to include time points in which the subject is standing still, but looking around. That said, we could still exclude those time points when estimating the place fields, since those are primarily relevant during navgiation, and not while standing still. 
  I do have that code, in the SpatialMapNew object. So, we could treat the differently, and estimate place fields from the SpatialMapNew code, and view fields (and conjunctions) from the JointMap.
  Also, in the literature, people do put some minimum speed thresholds on the spatial maps.
  This requires some rewrites, though, since currently I'm computing place and view maps using the same code, all using JointMap.
  
  OK, I think I managed to recovered the original place field for the cell p20181101s01a01g019c02. However, the conditioned view map looks quite different; I think I need to do the same same speed filtering when looking at the view map

  The one test that I have that the raytraced stuff actually work is if I plot the position for a trial obtained from `UnityData` and compare it to the same trial obtained using `UnityRaytraceData`, I get identical results (or at least close to identical). This means that there are no obvious errors in how I'm loading the raycast data. 
  My main issue now is that there are gaps in the raycast data that I don't know what to do with.
  The gaps appear to correspond to points going out of the screen

  By matching the actual eye trace between the eyelink data the raycast data, it looks like it is indeed aligned to trial onset, not cue onset. Points outside the screen are clipped, so at the very least we need to patch those. In addition, for the session of the cell I'm looking at, for some reason the raycast data is missing the first 95 ms of the first trial. 

  == 20260305

  Some updates. One issue that I had was that I assumed that the raycast data always started from the first point in trial. This is in generate not true. Another issue is that, if I use the actual gaze trace and match it to the trace from the eyelink file, the time point returned in the raycast file some times does not make sense. In particular, there appears to be some issue with trials in which the gaze starts out-of-bounds. 

  Something else is going on here. The raycast positions, as well as the actual spatial position, are all consistent; only the raw gaze is not. Is it possible that the raw gaez (i.e. in pixel coordinates on the screen), are somehow incorrect?


  I'm seeing shifts in the positions as well, actually. For instance, for trial number 356, there appears to be a 246ms shift between the raycast and the unity data. This doesn't seem to related to missing data; these incur a jump in time, which we can account for. The weird thing is that the first position lines up.
  
  Actually, I think we're good, at least to the point that I can reproduce almost exactly the place and view maps for the example cell I've been looking at.

  == 20260309

  We need to take the directionality of the place fields into account as well, that is is there a difference in firing in the fields depending on which way the animal is traversing them?

  The directional field analysis would need to access the `UnityData` object, and collect activity falling within the field, separated by the direction of travel. Should we define direction as the average direction traversed? It's possible that the direction changes within the field. We could even collect a histogram of travel directions falling within the field first. Find the point at which the animal entered the field, for each trial, and where it exited the field, and then just compute the vector connnecting the two points?

  == 20260310

  Thinking about how to incorporate the place field directionality in into the conjunction computation. The concern here is that we could get a different view-conditioned responses in the same place if the field itself is directional. For instance, we could get higher firing if we traverse a field from south to north, while looking right, compared to if we traverse from north to south, still looking to the right. If the cell did not care about where the eyes are looking, this assymmetry would still give rise to an apparent view field on the eastern part of the maze. How do we correct for that? If we condition on the traversal direction, then for a real gaze field, we should see the same regardless of the direction. 

  == 20260311

  #figure(image("../figures/p20181101s01_place_field_traversal_trial_16.png"),
          caption:[Example of incomplete field traversal. This is showing trial 16 for session p20181101s01, cell a01c019c02. One of this cell's 3 place fields is indicated in orange. The blue trajectory is the trajectory traversed by the subject in theis session, while the black dots represent the points along the trajectory at which spikes were by this cell.  The trajectory starts at the lower corner. For this particular example, the trajectory did not fully traverse the place field, but rather ended near its border. ]) <fig_place_field_traversal_example>


  @fig_place_field_traversal_example shows an example of a trial that was not included in the initial analysis, because the trajectory, represented by the blue line, did not fully traveser the field. In tihs version, I looked for the first point at which the trajectory entered the field, and the first point at which it exited the field, and the vector between these two points to compute the direction at which the traversal happened. I think this example should definitely be included, so perhaps I should just relax the criterion a bit? Rather than traversing the field fully, perhaps it is enough to traverse some portion of it, say half?


  That worked. Now most of the trials are used. Howeer, I'm still seeing major discrepancies between the view maps created simply by conditioning on the place map, and maps created by conditioning on the directionality as well. To be precise, I'm comparing the place conditioned map and the sum over the direction conditioned maps; these should be (nearly) identical, depending on how much data is left out by requiring full traversals of the fields. For instance, for the example considered in @fig_place_field_traversal_example, only a single spike was excluded when conditioning on direction. Thus, the view fields should be identical, but they are quite different.

  #figure(image("../figures/p20181101s01a01g019c02_place_conditioned_view_field.png"),
         caption: [View field conditioned on traversing the indicated place field from north to south (left) and from sout to north (right), as well as combined for both directions.] )


Another weird observation; if I sum up all the occupancy for the `JointOccupancy` map, including only bins that fall within the place field, I get a total of 115.59 seconds. However, if I sum up all the occupancy from the code that computes directional tuning, I get only 28.51 seconds. Where did all the rest go? Actually, this tracks with what the JointMap is giving, meaning both `JointMap` and `DirectionFiltered` contain the same number of spikes and the same occupancy. Something must be up with how the latter maps these to view coordinates.

OK, realized the mistake I was making. I was including bins for which there were spikes, even for the occupancy. I sort of fixed that, but now I have the opposite problem, i.e. that I have to long occupancy compared to the simple joint map. The hunt continues....
Perhaps what is going on is that we are essentially duplicating the value $l$ for the traversal direction across all the bins in the field. However, `jocc.weight` has already been aggrated across the entire trial. We actually cannot use `JointOccupancy` directly; we need to impelement a new type. We can probably still use `jocc.index`, though, since this basically just assigns each data point to a combination of view, place and head direction bin for each trial. In addition, this index has already been filtered for invalid bins, i.e. with low speed or excessive gaps.

 OK, fixed.


 == 20260312

 There is weird bug when I run `Hippocampus.issignificant` with `Hippocampus.SpatialInformationContent`. When I run it in batch mode, i.e. supplying a list of directories, I get signififance, but when I run it by first going into a specific directory, I do not get significance. It looks like two different objects are loaded, though I (appear to) supply the exact same arguments.
 Argh, it turns out I spelt "trial_start" differently 

 Now I realise I have computed all the fields starting from trial start, rather than cue onset. It looks like we get fewer significant cells if we exclude the cue period, so I'll re-run both spatial and gaze fields.
 

 I'm still not sure how exactly to correct for the directionality, though


 == 20260313

 Find cells with place fields near the center of the maze, since these fields would have higher potential for directionality.


 == 20260316

 14 cells with place field containing the center of the maze; `cellidx = [8, 13, 25, 38, 44, 68, 102, 103, 111, 119, 126, 135, 154, 165]`, of which 68 is the most promising, with a clear spatial field in the middle of the maze

I should create a plot for this with all cardinal directions.

Found another bug related to `trial_start`. Basically, since I wasn't using the appropriate `get_trial` method for the `UnityraytraceData` object in `DirectionalField`, I got an index mismatch, where the range `idx0:idx1`, referring to a trial starting at the cue, rather than navigation onset, extended beyond the length of the trial referring to `trial_start =2`. I've fixed it in the code, but I should make the constructor depend on `kwargs` like for the other objects.

Challenge: How to specify ranges that straddle zero? In order words, if I want to do a est-west traversal comparison, rather than north-south as I'm currently doing, I would need to look at direction from $-pi/2,pi/2$.  


== 20260317

Working on finding place fields with significant directionality. My current approach is to use a shuffle test, where I'm shuffling the relationship between estimated firing rate and directional bin. Some of the results look a bit fishy, so I'm looking into that. Does it make more sense to randomize at the level of single bins?
OK, randomizing at the single bin level give me something that looks a bit more reasonble. Running for all the spatially selective cells now.


Regarding conjunctions, it looks like 80 of 109 cells with both place and view fields come out as conjunctive. What does that look like in the contact of directionality of place fields? 126 out of 199 cells with place field are directional, 65 out of the 109 cells with both place and view selectivity have directional place fields.


This is still a lot of cells. I'm waiting for the results when I apply the restrictions on place and observations; hopefully I get fewer cells.



I also need to start thinking about that do talk about for the coffee and conversations presentation next week. 
I think I want to focus on the idea of the Hippocampus as encoding view points, and using those to navigate.
I'll add ome plots of mine, but also some of other people.

== 20260318

OK, I'm seeing more and more evidence that the place cells from the old Matlab code actually come from data aligned to the start of the trial and not the beginning of the go-cue. Of course, there could also be an error on my part here, but considering the amount of time I've spent goign through these things, I'm starting to think that maybe my approach is better. The issue I have with the Matlab code is that it is very dense and poorly documented. Objects like the `vmpv` do a lof of heavy lifting in large, monolithic functions, the logic of which is hard to follow. What I've tried to do is break up the logic into smaller steps, where each step can be handled separately. I'm hoping this will make things a bit more manageable, both now and in the future.

Something weird going on with the SIC code; the stored data doesn't seem to differentiate between the two different trial starts. This could because of changes I made a while back, so I think to be safe I'm going to recompute the spatial SIC for `trial start = 1.`

Once that completes, I can check whether the number of spatially selective cells is more similar to the one we had originally (though I suspect it won't be). At the very least, I can then check what causes the larger number of cells in the new version of the code.


#emph[Work on egocentric representation, i.e. looking left vs looking right conditioned on a particular traversal through a field].
This entails looking at the gaze coordinates in a say 60 degree wedge orthogonal to the direction of traversal. Thet is, if we are traversing a field from south to north, then we'll have one 60 degree wedge centered on east, representing gaze to the right, and another 60 degree wedge centered on west, representing gaze to the left. 

#canvas({
  line((-2.5, -2.5), (5, 0), ..arrow-style, name: "direction")
})


== 20260319

For egocentric gaze coding, I'll first get the egocentric gaze direction for each point a cell's place field for each traversal direction.


For some reason, the directionality of the fields of the model cells suddnely changed. In particular, the field between the top two pillars is now highly directional, whereas before it was not. I'm not sure which result is correct, but the previous plots actually look at bit more reasonably. I have to try and track down what changed.

== 20260323

For egocentric view coding, I have been using $60 degree$ cones centered on $plus.minus 90 degree$ with respect to the direction of traversal. 

#figure(image("../figures/p20180910s01a01g020c01_egocentric_view_coding.png"),
    caption: [Cell with a single place field between the west wall and the north-west pillar. As can be seen from the direcional plots in B, this cell is not directional, i.e. there is relationship between the firing rate of this cell within the place field and the direction at which the field is traversed for each trial. Nevertheless, according to D, the cell appears to only fire when the animal was looking left, with respect to the direction of traversal.])


I think we also need to worry about occupancy here. If the animal simply never looked right when traversing the place field shown in the above figure, we cannot really say that the cell prefers gazes to the left; those were simply all they were presented with. Thus, we should restrict the analysis to trials for which both left and rigth gazes were seen. Though, actually, if there were no rigth gazes, the $lambda_"right"$ vector would simply be empty. For this particular cell there were 137 bins representing the animal looking left, and 121 bins representing when the animal was looking right; it just so happens that all the 121 rigth bins had no spikes. 

#figure(image("../figures/p20180910s01a01g020c01_place_view_conjunction.png"),
          caption: [The conjunction between place and view for the same cell as in the previous figure. The fact that there is a conjunciton here is weird, because it is questionable whether the view field can even be seen from the place field. There seems to be an issue here where the outfield is smaller than the infield, mainly because the south part of the maze has less response, rather than the view-field having higher response. Should the test instead be whether the distribution of the infield response is more similar to the actual distribution in the view field, compared to a random patch that is not from the view field?])


#figure(image("../figures/p20180907s01a01g020c02_egocentric_view_coding.png"),
        caption: [Example of a cell with no difference between left and right views when traversing the place field.]))


After the meeting today, I question whether the way I'm computing tuning strength of the directional is correct. What I'm currently doing so to first assign a direction bin to each response based on the directionality of the current trial through the place field. Only place bins that actually fall within the place field are included. I then compute the circular mean vector of this, i.e.

$ mu_c = frac(sum_i lambda _i e^(j theta_i), sum_i lambda _i) $

where $lambda_i$ is the direction bin associated with firing rate $lambda_i$, and $j = sqrt(-1)$. To get a null distribution for this, I re-compute the firing rate $lambda$ for each directional bin by randomly associating each spike of the response with a directional value. It appears that I do this both for spike and for occupancy, actually. That means that the relationship between spiking and occupancy is maintained, and only the relationship between firing rate and direction is disrupted. Intuitively, the relation ship between occupancy and directionality is behavioural; does the subject spend more time going through the place field in one particular direction? This is separate from whether the spiking is sensitive to direction. So perhaps a better control is to only disrupt the relationship between spiking and directionality. Or we could just do the simple thing and disrupt the randomize the relatioship between mean firing rate and angle?

I was actually using a parametric approximation for the null distribution, where I assumed that the firing rates followed a Gamma distribution. That might not be entirely accurate.

Actually, neither fit very well, as the below figure shows

#figure(image("figures/gamma_vs_lognormal_fit.png"),
     caption: [The null distribution of tuning strength for one place field (shaded aread) compared to the fit to a Gamma distribution (blue line) and a Log-normal distribution (orange line). Neither of the two fit particularly well, and tend to overestimate the tails. Even a Beta distribution, which has restricted support between 0 and 1 does not really fit (grene line).])

== 20260324

Actually, I ended up using a kernel density estimator to estimate the p-values, which is exactly what the density plot above is using anyway. Now, both the blue and the orange place fields are directional, while the green field is not.

We may not need to compute directionality across the entire field, though, especially if we just want to define left and right. For that, we can find points during the traversal of the field where the local direction of movemement is along the preferred direction, and then define an (instantaneous) left/right directional dichotomy. 
What we are really interested in, though, is to what extent the observed conjunction between place and view can be attributed solely to place. That is, given that the place field is directional, if the animal is traversing the field in the preferred direction, that overtall firing rate will be higher than if it is traversing in the opposite direection. That means that whatever the animal is looking at during the preferred traversal also receives more spikes than what it is looking at when going in the opposite direction. If it is predominantly looking towards the view field, then this field could simply be a result of the spatial directionality, and not any interaction between place and view in itself. How do we check that? If the view field was a result of true interaction, then if we look at all the gaze points during the preferred traversal, we should see high firing for some of them, namely those in the view field, and low firing for others. If that's the case, the cell's firing pattern can not simply be explain by directional spatial preference. 

It looks like we have a total of 152 instances of being in the place field while traversing north, and 30 instances where the gaze point was in the identified view field. How many instances did we have when the gaze point was opposite to the view field? 
How should we do opposite? We could just draw a line from (each?) gaze point in the view field perpendicular to the local traversal direction, and then look at occupancy/spike count at the bin intersecting the maze on the other side of the traversal vector. Or we could do a cone; Basically, what I'm currently doing, but spefific to the view field. So, for spatial point in the field where the local traversal direction is parallel to the overall directionality of the field, identify points in the left and right view (cone) and then flag points where either the left or the right cone intersects with the view field. For those points, compare the firing rate on the cone intersecting the view field to the cone on the opposite side; if, acoss all flagged points, we see a significant difference in the two distribution of firing rates, then the view field cannot simply be a result of spatial effects. For alignment, let's just do a 60 degree cone centered on the overall directionality of the place field.

== 20260325

While working on the directionality analysis, I realised that the view field to the particular cell I'm focusing on, i.e. p20181101s01a01g019c02, is quite large. So, I wanted to see if I could concentrate it a bit by reducing the amount of smoothing and tightining the critiera for what constitutes a significant firing rate bin. Below is the resulting view field

#figure(image("../figures/p20181101s01a01g019c02_alternative_view_field.png"),
        caption: [Identified view field using 50 instead of 100 iterations for the Laplace smoothing and setting a threshold if p < 0.001 for accepting a bin as significantly above the surrogates.])


It turns out we still get significant conjunctions between the first two place fields and the view field, but not conjunction with the third place field, i.e. the one near the south-west pillar. What about the directionality?


#figure(image("../figures/p20181101s01a01g019c02_new_directional_analysis.png"),
       caption :[Directional analysis with refined view response fields. The orange dots, in particular, now represents the view field located in the North-East corner. The colormap represents the firing rate of the cell when the animal was traversing the place field (red outline) in the North-East direction.])


How many other cells like this do we see? But before that, how should we do the stats here? For this particular cell, we only have 4 bins with non-zero firing rate, all falling within the view field. 


Something interesting in #cite(<maoSpatialModulationHippocampal2021>), where they differentiate between spatial view (SV) and facing location (FL). What I'm currently looking at would be closer to SV, since we are also using actualy (x,y,z) coordinates do to our analysis. Facing location was, in their case, computed as the intersection between the head direction and the environment. They found more cells tuned to FL compared to SV (24% vs 11%).

In the same study, they note that they found evidence of eye-in-head coding, especially in the entorhinal cortex (EC). We actually haven't looked at this at all; this, I think would simply be coding in the 2D image space. For instance, are cells coding for eye movements towards particular parts of the visual field?

== 20260326

Crazy idea. Use Hebbian learning with Oja's rule to randomly explore the arean. Basically, we have a recurrent network of cells that receive input, in the form of gaze, and choose an action the space of move forward, turn left, turn right at each time step. If it reaches the reward zone set for that trial, it receives a reward signal that then feeds back to modulatory units (or all units?). At the beginning of the trial, the agent is "shown" an image of the poster to look for. Now, how do we want to represent this input? We could use some kind of camera. For instance, in we want to stick with Julia for this, we could use a normal 3D camera, as long as we set up the scene, we can get the releant projection matrices, and we can maniuplate the camera rotation and position. Camera movement = head direction. We do not really have a good model for eye poisition here. In general, at least for primates, moving the eyes around helps focus on details, while the overall picture represents lower resolution, overall views, or gist. 
I think this actually requires us to have some model of vision in place, which could get complicated. One simplifiction is to have the camera represent the head, and then have a filter that mimics the fovea applied to the camera image. Basically, the periphery would have quite low resolution, which only the central fixation giving details.


```julia
X = randn(40,40)
Z = Hippocampus.retinal_filter(40,40, 15,25)
fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,1])
heatmap!(ax1, X)
heatmap!(ax2, reshape(Z*X[:], 40,40))
```
#figure(image("../figures/retina_filter_illustration.png"),
       caption: [Illustration of how vision gets filtered through the retina. Left) original image, right image after being filtered with simulated retinal filter where the fovea is located at (15,25) in the visual field. The area around the fovea has high spatial resolution, while the periphery is only represented in a downsampled version.])


We would then have a larger portion of cells receiving input from the fovea than from the rest of the field. Since the fovea would move around the image (visual field), the mapping would also need to change. To simulate this, wwe could move the entire image, so that we can keep the mapping from filter pixels to cells constant. So, the camera represents the current visual field, for instance with a 120 degree viewing angle. To simulate eye movements, the retinal filter is then moved across the camera image. How do we deal with edge effects? We could avoid that by only allowing movement within the a central 60 degree window. In other words, the extreme edges of the camera image could never be foveated.


We could do it simpler, though. I always wanted to see what would happen if we have something simple, like some simple orgnaism that can only move forward or turn, and sense simple chemical gradients. We could have one sensor sensitive to a food-like chemical, i.e. if followed, it will be lead to food. Another would be noxious, i..e if followed, it would cause harm, and a third one could be neutral, i.e. nothing happens. Then, the model would need to learn how to follow the related gradient, stay away from the noxious gradient, and ignore the neutral gradient. All this would happen through Hebbian learning. In other words, we have a network that receives 3 inputs, correpsonding to the 3 different sensors, and ouputs 4 different actions; stay still, move forward, turn left or turn right. The environmental feedback would in the form of reward if the organism is close enough the source of the first type of gradient, and negative reward, i.e. punishment if it gets too close to the noxious stimulus. In addition, we could have physical feedback from the environment, for instance if the organism bumps into walls.


== 20260327

I computed the direcional tuning properties of each of the 198 cells that come out as spatially tuned. It turns out 127 of them (64%) show significant directional tunning for at least one of their spatial fields. 


== 20260104

After tightning the criteria from what constitutes place and view fields, I no longer see the strong view field in the upper right corner of the arena for the cell I've been looking at. There is still a strong view response there, but it only constitutes 2 bins, and so falls below the cluster threshold. The view field around the donkey poster persists, though. For this field, when conditioning on place and direction preferred direction (i.e. at an angle of around 270 degree ), there arent many (or any) views towars the view field. For this cell, that means that directional preference of the spatial field is not enough to explain the view response

#figure(image("../figures/presentation/p20181101s01a01g019c02_spatial_directional_view_field.png"),
        caption: [Left: View field for one cell. The heatmap shows the firing rate on each surface of the maze and the red dots indicate the part of the maze where the firing rate was significantly elevated. Right: The firing rate in each of the view bins, conditioned on the subject traversing the place field (orange dots) in the direction of the black arrow on the floor.]) <fig_spatial_directional_view_field_example>

There is still something a bit strange going on with the place and direction conditioned view fields. For example in right panel of @fig_spatial_directional_view_field_example, there is relatively high firing rate on the back of the north-west pillar, which is not visible from the place field.

Yup, there appears to be a bug in how direction conditioned response are computed. There are a bunch of place bin in `gidx.occupancy` that do not belong to the place field. Actually, that is not the explanation for what's going on. I might just have used different smoothing, resulting in a slightly wider field. The raw occupancy map is fine, though. Yup, it was just that I was using the smoothed, rather than the raw occupancy of set NaNs in the rate. It works just fine now.

For conjunctive coding, I think it makes more sense to compare a view map conditioned on random spatial parts compared to the view map conditioned on specfic place fields. For a a truly independent view map, it doesn't matter which part of space we condition on; the view map should look the same. One caveat here is that not all view locations are visible from al spatial locations, so when conditioning on random places, the resulting conditioned view responses should all be possible. In practice, that means to sample locations that occur together view a view location.

We can do this by repeatedly sampling from the non-covered locations and computing the firing rate in the view field.

For the above example, I do, however find that when conditioning on other parts of the floor, I get much higher performance within the view field than if I condition on an actual place field.

Might we not use joint SIC to ascertain whether a cell is conjunctive?

Basically, what we are currently doing is computing the mutual information between spiking activity and space. Intuitively, we want to know what affects the spiking activity of the cell. One possibility is that spiking activity is solely influenced by the time spent per bin; that is, more time spent, more spikes.

$ I(X;S) = H(X) + H(S) - H(X,S) $

We also have occupancy on the view bins. What if there is no information about space? If he subject spends an equal amount of time in each bin, then the occupancy itself does not convey any information about space. If the spiking activity is not uniform, then, trivially, spiking activity conveys more information about space than does occupancy alone. If we assume that each bin in itself is equally likely, then we can write

$ P(X,S) = sum_i P(X,S|b=i) P(b=i) $

If $lambda_j$ is the firing rate in bin $j$, then the probability of having a spike in that bin given an occupancy of $delta t$ is

$ P(S=1 | X=j) = lambda_j delta t $

The probabililty of a spike regardless of bin is further given by

$ P(S=1) = lambda delta t $

where $lambda$ is the overall firing rate, 
given by 

$ lambda = sum_i lambda_i p_i $,

where $p_i$ is the occupancy probability of bin $i$.

The mutual information between spiking activity S and binned space X. 

If we have both space bins $X$ and view bins $Y$, the probability we'd be interested in is

$ P(S|X,Y) $

Analogous to the above, we can now write

$ P(S=1|x=i, y=j) = lambda_(i j) delta t $

where $lambda_(i j)$ is the firing rate when the subject was in spatial bin $i$ looking at view bin $j$. Let's look at the difference

$ H(S|X) + H(S|Y) - H(S|X,Y) $

$ H(X|S) + H(Y|S) - H(X,Y|S) = I(X;Y|S) $

$ I(S|X;Y) = sum_(s,i,j) P(S=s|X=i,Y=j) P(X=i,Y=j) log frac(P(S=s|X=i,Y=j),P(X=i,Y=j)) $ 

== 20260402

Ideally, I would like an expression that just uses the SIC we have already computed, i.e. the difference between the SIC for the joint minus the information for the SIC for each of the other variables. Ideally, again, we want a measue that is zero if the mutual information between spiking activity and place and view is the same as the sum of mutual information between spiking activity and space and spiking activity and view.

Since conditioning always reduces entropy, we have

$ H(S|X,Y) <= H(S|X) + H(S|Y) $

== 20260408

Just thikning about what the next steps are. I think the approach I want to take is to use SIC to determine whether a cell is spatially selective, then look at the number of response fields for each cell; a cell that is spatially selective, but has no spatial fields, would not be a place cell. Then, we look at the cells with at least one place field, and we find the proportion of those cells that also have significant view responses. We then determine the proportion of view responses that can be explained by just spatial response fields with a certain direcional preference. Among the remaining cells, i.e. those that have spatial response fields, and gaze response fields that cannot be explained solely by the spatial responses, we look for cells that are conjunctive. How do we determine whethre a cell is conjunctive? I think the most straight forward way is to compare view responses conditioned on the place field, to view resposnes conditioned on a part of space that is not part of the place field. For a conjunctive cell, the view responses should be stronger when conditioning on the actual place field. This, however, is not the only way that a cell can be conjunctive. We could also have a more general interaction between place and view. We can quantify that by looking at how much mutual information there is about place and view jointly, and the spiking activity, vs how much mutural information there is between spikes and place, and spikes and view, separately, i.e.

$ "Syn"(P,G;S) = I(P,G;S) - (I(P:S) + I(G;S)) $

If $"Syn"(P,G:S)$ is positive, that indicates a synergistc relationship between spiking, place and view.

#figure(image("../figures/p20181101s01a01g019c02_place_response_comparison.png"),
        caption: [Comparison of place response for two sets of paramters. A) Place response after filtering bins only based on speed. B) In adition, include only bins where the total place response was 50ms, the total view response was 20ms, and there were at least 5 trials where the bin exceeded those durations.]) <fig_spatial_response_comparison>

Looking at @fig_spatial_response_comparison, it is a bit curious that filtering out bins with transient place/view responses actually increases the occupancy (top left panel in A vs B). It would make more sense if the total duration went down. Actually, that was bug. I was using the old SpatialMap, which didn't properly account for both place and view bins. After using `JointMap` to compute the spatial map, the results are more reasonable. 

The main effect appears to be slightly reducing the size of the place field. What about the view responses? 

I was tracking what appeared to be a bug, but I think I've concluded it was just an optical illusion. It looked like a view bin (bin number 526, on the north wall of the north-east pillar) had a higher spike count after conditioning, but actually the spike counts were identical (37).

Back to the pipeline I mentioned above, let's figure out the number of cells with at least one place field and at least one view field


#strike[I just noticed that if I set the number of iterations in the smoothing to 50, I only get 2 cells where the SIC for view is significant, whereas if I use 100 iterations, I get 122 cells with significant view responses. I'm trying to use 75 iterations to see what effect that his. This was prompted by the observation that the 4th view selective cell when using 100 iterations is not selective when using 50 iterations, even though I can identify two distinct view fields.]#emph[This was incorrect; I had simply not computed the SIc for 50 iterations, except for 2 cells. I'm currently computing those]. 
That raises a question; does a cell not be overall spatially selective (i.e. have an SIC above the surrogates) in order to be considered a place cell?


== 20260409

We have spatially selective cells, view selective cells, pure place cells, pure view cells (?), place and view cells, conjunctive cells.
Come up with a schematic for showing how these are related. 

== 20260410

I find that more than 90% of the cells with spatial fields have at least one field that is directional. Does that mean that I'm too liberal in what I define as a directional field?

== 20260413

After today's meeting, we decided to try and be more restrice when determining directioality, and only look at the cardinal directions, rather than 24 directions as I'm currently doing. We could just do that for the data we currently have; Instead of 24 locations, just group into $60 deg$ cones around the four canonical directions. What about statistics? The circular mean is not that useful when we only have 4 directions. We just want to know whether the firing associated with a particlar direction is higher. Can we sub-sample for statistics?

== 20260414

I needed to rewrite the `DirectionFiltered` type a bit because it did not retain information about number of spikes and occupancy on a per trial basis. I have have a function `get_cardinal_direction_tuning` that takes in an index, an occupancy vector and a vector of spike counts and computes the firing rate, per trial, for each of the cardinal directions. How do we establish significance? I think shuffling is again the best approach. We can assign a random spike count, sampled from the spike count vector, to a particular cardinal direction bin. The question is then whether the actual distributino of firing rates across the 4 cardinal bins is different from the distributions resulting from randomly assigning spike counts. 
 #bibliography("bibliography.bib")