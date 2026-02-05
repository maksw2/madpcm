@echo off

:: --- MADPCM SLOW ---
test.exe encode slow Egyptischer_Marsch.wav Egyptischer_Marsch_ad.wav
test.exe encode slow Marche_Persanne.wav Marche_Persanne_ad.wav
test.exe encode slow suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
test.exe encode slow the_four_seasons.wav the_four_seasons_ad.wav

:: Decode
test.exe decode Egyptischer_Marsch_ad.wav Egyptischer_Marsch_d.wav
test.exe decode Marche_Persanne_ad.wav Marche_Persanne_d.wav
test.exe decode suppe_poet_and_peasant_ad.wav suppe_poet_and_peasant_d.wav
test.exe decode the_four_seasons_ad.wav the_four_seasons_d.wav

:: Verify Slow
echo --- MADPCM SLOW RESULTS ---
python ver.py Egyptischer_Marsch.wav Egyptischer_Marsch_d.wav
python ver.py Marche_Persanne.wav Marche_Persanne_d.wav
python ver.py suppe_poet_and_peasant.wav suppe_poet_and_peasant_d.wav
python ver.py the_four_seasons.wav the_four_seasons_d.wav

:: --- MADPCM FAST ---
test.exe encode fast Egyptischer_Marsch.wav Egyptischer_Marsch_ad.wav
test.exe encode fast Marche_Persanne.wav Marche_Persanne_ad.wav
test.exe encode fast suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
test.exe encode fast the_four_seasons.wav the_four_seasons_ad.wav

:: Decode
test.exe decode Egyptischer_Marsch_ad.wav Egyptischer_Marsch_d.wav
test.exe decode Marche_Persanne_ad.wav Marche_Persanne_d.wav
test.exe decode suppe_poet_and_peasant_ad.wav suppe_poet_and_peasant_d.wav
test.exe decode the_four_seasons_ad.wav the_four_seasons_d.wav

:: Verify Fast
echo --- MADPCM FAST RESULTS ---
python ver.py Egyptischer_Marsch.wav Egyptischer_Marsch_d.wav
python ver.py Marche_Persanne.wav Marche_Persanne_d.wav
python ver.py suppe_poet_and_peasant.wav suppe_poet_and_peasant_d.wav
python ver.py the_four_seasons.wav the_four_seasons_d.wav