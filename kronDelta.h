void kronDeltaFunc(int &num, int &countNum, int &condNum)
{
	if ( countNum != condNum)
	{
		num = 0;
	}
	else if ( countNum == condNum)
	{
		num = 1;
	}
	else
	{
		cout << "error in kronDeltaFunc" << endl;
	}
	return;
}
