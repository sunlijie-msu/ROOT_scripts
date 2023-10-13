// {
// 	char cmd[128];
// 	for (int i = 1; i <= 131; i++) {
// 		sprintf(cmd, "hadd /mnt/analysis/e21010/sorted/run_%03d.root /mnt/analysis/e21010/sorted/%03d_*.root", i, i);
// 		gSystem->Exec(cmd);
// 	}
// }

{
	char cmd[128];
	for (int i = 1; i <= 131; i++) {
		sprintf(cmd, "hadd /mnt/analysis/e21010/unpacked/raw_%03d.root /mnt/analysis/e21010/unpacked/%03d_*.root", i, i);
		gSystem->Exec(cmd);
	}
}
