void kronDeltaConstFunc(int &num, double &countNum, double &condNum)
{
	if ( countNum < condNum)
	{
		num = 0;
	}
	else if ( countNum >= condNum)
	{
		num = 1;
	}
	else
	{
		cout << "error in kronDeltaConstFunc" << endl;
	}
	return;
}
