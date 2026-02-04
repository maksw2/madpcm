@echo off

:: --- IMA ADPCM (ffmpeg) ---
ffmpeg -y -hide_banner -i Egyptischer_Marsch.wav -c:a adpcm_ima_wav Egyptischer_Marsch_ima.wav
ffmpeg -y -hide_banner -i Marche_Persanne.wav -c:a adpcm_ima_wav Marche_Persanne_ima.wav
ffmpeg -y -hide_banner -i suppe_poet_and_peasant.wav -c:a adpcm_ima_wav suppe_poet_and_peasant_ima.wav
ffmpeg -y -hide_banner -i the_four_seasons.wav -c:a adpcm_ima_wav the_four_seasons_ima.wav

ffmpeg -y -hide_banner -i Egyptischer_Marsch_ima.wav Egyptischer_Marsch_ima_d.wav
ffmpeg -y -hide_banner -i Marche_Persanne_ima.wav Marche_Persanne_ima_d.wav
ffmpeg -y -hide_banner -i suppe_poet_and_peasant_ima.wav suppe_poet_and_peasant_ima_d.wav
ffmpeg -y -hide_banner -i the_four_seasons_ima.wav the_four_seasons_ima_d.wav

:: --- ADPCM-XQ ---
adpcm-xq.exe -8 -y Egyptischer_Marsch.wav Egyptischer_Marsch_xq.wav
adpcm-xq.exe -8 -y Marche_Persanne.wav Marche_Persanne_xq.wav
adpcm-xq.exe -8 -y suppe_poet_and_peasant.wav suppe_poet_and_peasant_xq.wav
::adpcm-xq.exe -8 -y the_four_seasons.wav the_four_seasons_xq.wav

adpcm-xq.exe -d -y Egyptischer_Marsch_xq.wav Egyptischer_Marsch_xq_d.wav
adpcm-xq.exe -d -y Marche_Persanne_xq.wav Marche_Persanne_xq_d.wav
adpcm-xq.exe -d -y suppe_poet_and_peasant_xq.wav suppe_poet_and_peasant_xq_d.wav
::adpcm-xq.exe -d -y the_four_seasons_xq.wav the_four_seasons_xq_d.wav

:: --- MADPCM SLOW ---
bin\x64\Release\example.exe encode slow Egyptischer_Marsch.wav Egyptischer_Marsch_ad.wav
bin\x64\Release\example.exe encode slow Marche_Persanne.wav Marche_Persanne_ad.wav
bin\x64\Release\example.exe encode slow suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
bin\x64\Release\example.exe encode slow the_four_seasons.wav the_four_seasons_ad.wav

:: Decode
bin\x64\Release\example.exe decode Egyptischer_Marsch_ad.wav Egyptischer_Marsch_d.wav
bin\x64\Release\example.exe decode Marche_Persanne_ad.wav Marche_Persanne_d.wav
bin\x64\Release\example.exe decode suppe_poet_and_peasant_ad.wav suppe_poet_and_peasant_d.wav
bin\x64\Release\example.exe decode the_four_seasons_ad.wav the_four_seasons_d.wav

:: Verify Slow
echo --- MADPCM SLOW RESULTS ---
python ver.py Egyptischer_Marsch.wav Egyptischer_Marsch_d.wav
python ver.py Marche_Persanne.wav Marche_Persanne_d.wav
python ver.py suppe_poet_and_peasant.wav suppe_poet_and_peasant_d.wav
python ver.py the_four_seasons.wav the_four_seasons_d.wav

:: --- MADPCM FAST ---
bin\x64\Release\example.exe encode fast Egyptischer_Marsch.wav Egyptischer_Marsch_ad.wav
bin\x64\Release\example.exe encode fast Marche_Persanne.wav Marche_Persanne_ad.wav
bin\x64\Release\example.exe encode fast suppe_poet_and_peasant.wav suppe_poet_and_peasant_ad.wav
bin\x64\Release\example.exe encode fast the_four_seasons.wav the_four_seasons_ad.wav

:: Decode
bin\x64\Release\example.exe decode Egyptischer_Marsch_ad.wav Egyptischer_Marsch_d.wav
bin\x64\Release\example.exe decode Marche_Persanne_ad.wav Marche_Persanne_d.wav
bin\x64\Release\example.exe decode suppe_poet_and_peasant_ad.wav suppe_poet_and_peasant_d.wav
bin\x64\Release\example.exe decode the_four_seasons_ad.wav the_four_seasons_d.wav

:: Verify Fast
echo --- MADPCM FAST RESULTS ---
python ver.py Egyptischer_Marsch.wav Egyptischer_Marsch_d.wav
python ver.py Marche_Persanne.wav Marche_Persanne_d.wav
python ver.py suppe_poet_and_peasant.wav suppe_poet_and_peasant_d.wav
python ver.py the_four_seasons.wav the_four_seasons_d.wav