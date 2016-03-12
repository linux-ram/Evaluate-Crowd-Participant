# Evaluate-Crowd-Participant
The code presents the novel approach to identify unhelpful or adversarial participants in a crowdsourcing setting where multiple labelers are available. Further, the code demonstrates the robustness of the perfected scoring scheme to evaluate the participants.<br><br>
Note:<br><br>
(a) Constant features `beta` and `gamma` are the last dimension of `alpha` and `W` respectively.
(b) Resulting probabilities [0, 1] from the model shall be mapped to [-1 1] using the simple transformation,  2 **x** y{0,1} - 1 = y{-1,1}.<br><br>
Folders `ErrorStatistics`, `MultiLabelerMethods` and the `DATA` were shared by Yan Yan, now a Senior Software Engineer at LinkedIn and Romer Rosales, Principal Scientist at LinkedIn.<br><br>
Please feel free to use the code in your research and development works. We would appreciate a citation to the paper below when this code is helpful in obtaining results in your future publications.<br><br>
**Publication for citation:**<br>
(1) R. Subramanian, Romer Rosales, Glenn Fung, Prof. Jennifer Dy, "Evaluating Crowdsourcing Participants in the Absence of Ground-Truth", at the 26th Annual Conference on Neural Information Processing Systems (NIPS 2012) Workshop on Human Computation for Science and Computational Sustainability,  Accepted Date: Oct 7, 2012.<br><br>
(2) Yan Yan, Rómer Rosales, Glenn Fung, Ramanathan Subramanian, and Jennifer Dy, "Learning from Multiple Annotators with Varying Expertise", Machine Learning Journal, June 2014, Volume 95, Issue 3, pp 291-327.
