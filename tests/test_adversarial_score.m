function tests = test_adversarial_score
% Numeric check of adversarial.m (equation 19).
% From the repo root: matlab -batch "runtests('tests')"
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename("fullpath")));
addpath(root);
testCase.TestData.root = root;
end

function test_known_scores(testCase)
Pz_x_b = 0.25;
Py_xzbar = [0.9 0.2];
Py_xz = [0.1 0.8];
[p_ann, score] = adversarial(Pz_x_b, Py_xzbar, Py_xz, 2);
verifyEqual(testCase, p_ann, [0.4428571428571429 0.22142857142857145], "AbsTol", 1e-12);
verifyEqual(testCase, score, [0.8145080375642126 1.507655218124158], "AbsTol", 1e-12);
end

function test_uniform_labels_score_zero(testCase)
[p_ann, score] = adversarial(0.5, [1 1], [1 1], 2);
verifyEqual(testCase, p_ann, [1 1], "AbsTol", 1e-12);
verifyEqual(testCase, score, [0 0], "AbsTol", 1e-12);
end
