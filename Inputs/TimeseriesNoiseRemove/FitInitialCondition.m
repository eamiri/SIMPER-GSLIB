clear
clc

x = [0;0.0500000000000000;0.100000000000000;0.150000000000000;0.200000000000000;0.250000000000000;0.300000000000000;0.350000000000000;0.450000000000000;0.550000000000000;0.650000000000000;1];
% y = [8;9.30000000000000;7.10000000000000;6.30000000000000;6;5.70000000000000;5.20000000000000;4.50000000000000;3.70000000000000;2.80000000000000;2];
y = [8.27110523150057,9.86011367908354,7.22913520325939,6.32721298168861,5.99579880842830,5.67170756112393,5.17181915522090,4.43102574201201,3.61367021657010,2.72360604075186,1.91454662316467,-2]';


[fitresult, gof] = createFit(x,y);

NDATA=[0
0.002
0.004
0
0.002
0.006
0.004
0.006
0.008
0.008
0
0.01
0.002
0.004
0.01
0.006
0.012
0.008
0.012
0.014
0.01
0.014
0.012
0.016
0.016
0.014
0.018
0.018
0.016
0.02
0.018
0.02
0.022
0.02
0.022
0.024
0.022
0.024
0.024
0.026
0.026
0.026
0.028
0.028
0.028
0.03
0.03
0.03
0.032
0.032
0.032
0.034
0.034
0.034
0.036
0.036
0.036
0.038
0.038
0.038
0.04
0.04
0.04
0.042
0.042
0.042
0.044
0.044
0.044
0.046
0.046
0.046
0.048
0.048
0.048
0.05
0.05
0.05
0.052
0.052
0.052
0.054
0.054
0.054
0.056
0.056
0.056
0.058
0.058
0.058
0.06
0.06
0.06
0.062
0.062
0.062
0.064
0.064
0.064
0.066
0.066
0.066
0.068
0.068
0.068
0.07
0.07
0.07
0.072
0.072
0.072
0.074
0.074
0.074
0.076
0.076
0.076
0.078
0.078
0.078
0.08
0.08
0.08
0.082
0.082
0.082
0.084
0.084
0.084
0.086
0.086
0.086
0.088
0.088
0.088
0.09
0.09
0.09
0.092
0.092
0.092
0.094
0.094
0.094
0.096
0.096
0.096
0.098
0.098
0.098
0.1
0.1
0.1
0.102
0.102
0.102
0.104
0.104
0.104
0.106
0.106
0.106
0.108
0.108
0.108
0.11
0.11
0.11
0.112
0.112
0.112
0.114
0.114
0.114
0.116
0.116
0.116
0.118
0.118
0.118
0.12
0.12
0.12
0.122
0.122
0.122
0.124
0.124
0.124
0.126
0.126
0.126
0.128
0.128
0.128
0.13
0.13
0.13
0.132
0.132
0.132
0.134
0.134
0.134
0.136
0.136
0.136
0.138
0.138
0.138
0.14
0.14
0.14
0.142
0.142
0.142
0.144
0.144
0.144
0.146
0.146
0.146
0.148
0.148
0.148
0.15
0.15
0.15
0.152
0.152
0.152
0.154
0.154
0.154
0.156
0.156
0.156
0.158
0.158
0.158
0.16
0.16
0.16
0.162
0.162
0.162
0.164
0.164
0.164
0.166
0.166
0.166
0.168
0.168
0.168
0.17
0.17
0.17
0.172
0.172
0.172
0.174
0.174
0.174
0.176
0.176
0.176
0.178
0.178
0.178
0.18
0.18
0.18
0.182
0.182
0.182
0.184
0.184
0.184
0.186
0.186
0.186
0.188
0.188
0.188
0.19
0.19
0.19
0.192
0.192
0.192
0.194
0.194
0.194
0.196
0.196
0.196
0.198
0.198
0.198
0.2
0.2
0.2
0.202
0.202
0.202
0.204
0.204
0.204
0.206
0.206
0.206
0.208
0.208
0.208
0.21
0.21
0.21
0.212
0.212
0.212
0.214
0.214
0.214
0.216
0.216
0.216
0.218
0.218
0.218
0.22
0.22
0.22
0.222
0.222
0.222
0.224
0.224
0.224
0.226
0.226
0.226
0.228
0.228
0.228
0.23
0.23
0.23
0.232
0.232
0.232
0.234
0.234
0.234
0.236
0.236
0.236
0.238
0.238
0.238
0.24
0.24
0.24
0.242
0.242
0.242
0.244
0.244
0.244
0.246
0.246
0.246
0.248
0.248
0.248
0.25
0.25
0.25
0.252
0.252
0.252
0.254
0.254
0.254
0.256
0.256
0.256
0.258
0.258
0.258
0.26
0.26
0.26
0.262
0.262
0.262
0.264
0.264
0.264
0.266
0.266
0.266
0.268
0.268
0.268
0.27
0.27
0.27
0.272
0.272
0.272
0.274
0.274
0.274
0.276
0.276
0.276
0.278
0.278
0.278
0.28
0.28
0.28
0.282
0.282
0.282
0.284
0.284
0.284
0.286
0.286
0.286
0.288
0.288
0.288
0.29
0.29
0.29
0.292
0.292
0.292
0.294
0.294
0.294
0.296
0.296
0.296
0.298
0.298
0.298
0.3
0.3
0.3
0.302
0.302
0.302
0.304
0.304
0.304
0.306
0.306
0.306
0.308
0.308
0.308
0.31
0.31
0.31
0.312
0.312
0.312
0.314
0.314
0.314
0.316
0.316
0.316
0.318
0.318
0.318
0.32
0.32
0.32
0.322
0.322
0.322
0.324
0.324
0.324
0.326
0.326
0.326
0.328
0.328
0.328
0.33
0.33
0.33
0.332
0.332
0.332
0.334
0.334
0.334
0.336
0.336
0.336
0.338
0.338
0.338
0.34
0.34
0.34
0.342
0.342
0.342
0.344
0.344
0.344
0.346
0.346
0.346
0.348
0.348
0.348
0.35
0.35
0.35
0.352
0.352
0.352
0.354
0.354
0.354
0.356
0.356
0.356
0.358
0.358
0.358
0.36
0.36
0.36
0.362
0.362
0.362
0.364
0.364
0.364
0.366
0.366
0.366
0.368
0.368
0.368
0.37
0.37
0.37
0.372
0.372
0.372
0.374
0.374
0.374
0.376
0.376
0.376
0.378
0.378
0.378
0.38
0.38
0.38
0.382
0.382
0.382
0.384
0.384
0.384
0.386
0.386
0.386
0.388
0.388
0.388
0.39
0.39
0.39
0.392
0.392
0.392
0.394
0.394
0.394
0.396
0.396
0.396
0.398
0.398
0.398
0.4
0.4
0.4
0.402
0.402
0.402
0.404
0.404
0.404
0.406
0.406
0.406
0.408
0.408
0.408
0.41
0.41
0.41
0.412
0.412
0.412
0.414
0.414
0.414
0.416
0.416
0.416
0.418
0.418
0.418
0.42
0.42
0.42
0.422
0.422
0.422
0.424
0.424
0.424
0.426
0.426
0.426
0.428
0.428
0.428
0.43
0.43
0.43
0.432
0.432
0.432
0.434
0.434
0.434
0.436
0.436
0.436
0.438
0.438
0.438
0.44
0.44
0.44
0.442
0.442
0.442
0.444
0.444
0.444
0.446
0.446
0.446
0.448
0.448
0.448
0.45
0.45
0.45
0.452
0.452
0.452
0.454
0.454
0.454
0.456
0.456
0.456
0.458
0.458
0.458
0.46
0.46
0.46
0.462
0.462
0.462
0.464
0.464
0.464
0.466
0.466
0.466
0.468
0.468
0.468
0.47
0.47
0.47
0.472
0.472
0.472
0.474
0.474
0.474
0.476
0.476
0.476
0.478
0.478
0.478
0.48
0.48
0.48
0.482
0.482
0.482
0.484
0.484
0.484
0.486
0.486
0.486
0.488
0.488
0.488
0.49
0.49
0.49
0.492
0.492
0.492
0.494
0.494
0.494
0.496
0.496
0.496
0.498
0.498
0.498
0.5
0.5
0.5
0.502
0.502
0.502
0.504
0.504
0.504
0.506
0.506
0.506
0.508
0.508
0.508
0.51
0.51
0.51
0.512
0.512
0.512
0.514
0.514
0.514
0.516
0.516
0.516
0.518
0.518
0.518
0.52
0.52
0.52
0.522
0.522
0.522
0.524
0.524
0.524
0.526
0.526
0.526
0.528
0.528
0.528
0.53
0.53
0.53
0.532
0.532
0.532
0.534
0.534
0.534
0.536
0.536
0.536
0.538
0.538
0.538
0.54
0.54
0.54
0.542
0.542
0.542
0.544
0.544
0.544
0.546
0.546
0.546
0.548
0.548
0.548
0.55
0.55
0.55
0.552
0.552
0.552
0.554
0.554
0.554
0.556
0.556
0.556
0.558
0.558
0.558
0.56
0.56
0.56
0.562
0.562
0.562
0.564
0.564
0.564
0.566
0.566
0.566
0.568
0.568
0.568
0.57
0.57
0.57
0.572
0.572
0.572
0.574
0.574
0.574
0.576
0.576
0.576
0.578
0.578
0.578
0.58
0.58
0.58
0.582
0.582
0.582
0.584
0.584
0.584
0.586
0.586
0.586
0.588
0.588
0.588
0.59
0.59
0.59
0.592
0.592
0.592
0.594
0.594
0.594
0.596
0.596
0.596
0.598
0.598
0.598
0.6
0.6
0.6
0.602
0.602
0.602
0.604
0.604
0.604
0.606
0.606
0.606
0.608
0.608
0.608
0.61
0.61
0.61
0.612
0.612
0.612
0.614
0.614
0.614
0.616
0.616
0.616
0.618
0.618
0.618
0.62
0.62
0.62
0.622
0.622
0.622
0.624
0.624
0.624
0.626
0.626
0.626
0.628
0.628
0.628
0.63
0.63
0.63
0.632
0.632
0.632
0.634
0.634
0.634
0.636
0.636
0.636
0.638
0.638
0.638
0.64
0.64
0.64
0.642
0.642
0.642
0.644
0.644
0.644
0.646
0.646
0.646
0.648
0.648
0.648
0.65
0.65
0.65
0.652
0.652
0.652
0.654
0.654
0.654
0.656
0.656
0.656
0.658
0.658
0.658
0.66
0.66
0.66
0.662
0.662
0.662
0.664
0.664
0.664
0.666
0.666
0.666
0.668
0.668
0.668
0.67
0.67
0.67
0.672
0.672
0.672
0.674
0.674
0.674
0.676
0.676
0.676
0.678
0.678
0.678
0.68
0.68
0.68
0.682
0.682
0.682
0.684
0.684
0.684
0.686
0.686
0.686
0.688
0.688
0.688
0.69
0.69
0.69
0.692
0.692
0.692
0.694
0.694
0.694
0.696
0.696
0.696
0.698
0.698
0.698
0.7
0.7
0.7
0.702
0.702
0.702
0.704
0.704
0.704
0.706
0.706
0.706
0.708
0.708
0.708
0.71
0.71
0.71
0.712
0.712
0.712
0.714
0.714
0.714
0.716
0.716
0.716
0.718
0.718
0.718
0.72
0.72
0.72
0.722
0.722
0.722
0.724
0.724
0.724
0.726
0.726
0.726
0.728
0.728
0.728
0.73
0.73
0.73
0.732
0.732
0.732
0.734
0.734
0.734
0.736
0.736
0.736
0.738
0.738
0.738
0.74
0.74
0.74
0.742
0.742
0.742
0.744
0.744
0.744
0.746
0.746
0.746
0.748
0.748
0.748
0.75
0.75
0.75
0.752
0.752
0.752
0.754
0.754
0.754
0.756
0.756
0.756
0.758
0.758
0.758
0.76
0.76
0.76
0.762
0.762
0.762
0.764
0.764
0.764
0.766
0.766
0.766
0.768
0.768
0.768
0.77
0.77
0.77
0.772
0.772
0.772
0.774
0.774
0.774
0.776
0.776
0.776
0.778
0.778
0.778
0.78
0.78
0.78
0.782
0.782
0.782
0.784
0.784
0.784
0.786
0.786
0.786
0.788
0.788
0.788
0.79
0.79
0.79
0.792
0.792
0.792
0.794
0.794
0.794
0.796
0.796
0.796
0.798
0.798
0.798
0.8
0.8
0.8
0.802
0.802
0.802
0.804
0.804
0.804
0.806
0.806
0.806
0.808
0.808
0.808
0.81
0.81
0.81
0.812
0.812
0.812
0.814
0.814
0.814
0.816
0.816
0.816
0.818
0.818
0.818
0.82
0.82
0.82
0.822
0.822
0.822
0.824
0.824
0.824
0.826
0.826
0.826
0.828
0.828
0.828
0.83
0.83
0.83
0.832
0.832
0.832
0.834
0.834
0.834
0.836
0.836
0.836
0.838
0.838
0.838
0.84
0.84
0.84
0.842
0.842
0.842
0.844
0.844
0.844
0.846
0.846
0.846
0.848
0.848
0.848
0.85
0.85
0.85
0.852
0.852
0.852
0.854
0.854
0.854
0.856
0.856
0.856
0.858
0.858
0.858
0.86
0.86
0.86
0.862
0.862
0.862
0.864
0.864
0.864
0.866
0.866
0.866
0.868
0.868
0.868
0.87
0.87
0.87
0.872
0.872
0.872
0.874
0.874
0.874
0.876
0.876
0.876
0.878
0.878
0.878
0.88
0.88
0.88
0.882
0.882
0.882
0.884
0.884
0.884
0.886
0.886
0.886
0.888
0.888
0.888
0.89
0.89
0.89
0.892
0.892
0.892
0.894
0.894
0.894
0.896
0.896
0.896
0.898
0.898
0.898
0.9
0.9
0.9
0.902
0.902
0.902
0.904
0.904
0.904
0.906
0.906
0.906
0.908
0.908
0.908
0.91
0.91
0.91
0.912
0.912
0.912
0.914
0.914
0.914
0.916
0.916
0.916
0.918
0.918
0.918
0.92
0.92
0.92
0.922
0.922
0.922
0.924
0.924
0.924
0.926
0.926
0.926
0.928
0.928
0.928
0.93
0.93
0.93
0.932
0.932
0.932
0.934
0.934
0.934
0.936
0.936
0.936
0.938
0.938
0.938
0.94
0.94
0.94
0.942
0.942
0.942
0.944
0.944
0.944
0.946
0.946
0.946
0.948
0.948
0.948
0.95
0.95
0.95
0.952
0.952
0.952
0.954
0.954
0.954
0.956
0.956
0.956
0.958
0.958
0.958
0.96
0.96
0.96
0.962
0.962
0.962
0.964
0.964
0.964
0.966
0.966
0.966
0.968
0.968
0.968
0.97
0.97
0.97
0.972
0.972
0.972
0.974
0.974
0.974
0.976
0.976
0.976
0.978
0.978
0.978
0.98
0.98
0.98
0.982
0.982
0.982
0.984
0.984
0.984
0.986
0.986
0.986
0.988
0.988
0.988
0.99
0.99
0.99
0.992
0.992
0.992
0.994
0.994
0.994
0.996
0.996
0.996
0.998
0.998
0.998
1
1
1];
fitIC=zeros(size(NDATA,1),1);
for i=1:size(NDATA,1)
    fitIC(i,1)=fitresult(1-NDATA(i,1));
end