void kronDeltaPulseFunc(int &num, int &countNum, int &condNum)
{
	if ( countNum < condNum)
	{
		num = 0;
	}
	else if ( countNum == condNum)
	{
		num = 1;
	}
	else if ( countNum > condNum)
	{
		num = 0;
	}
	else
	{
		cout << "error in kronDeltaPulseFunc" << endl;
	}
	return;
}
