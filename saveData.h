void saveData(string name, MatrixXd matrix)
{
        const static IOFormat CSVFormat(FullPrecision, DontAlignCols, ",", "\n");
        ofstream file(name);
        if(file.is_open())
        {
                file << matrix.format(CSVFormat);
                file.close();
        }
}
