with open("/home/yhr8/Documents/MLVA_Search/VNTR_chart_edit_4.txt", 'r') as fin:
    with open("/home/yhr8/Documents/MLVA_Search/VNTR_chart_edit_5.txt", 'w') as fout:
        '''
        for line in fin:
            fout.write(line.replace(',', ':[', 1))

        for line in fin:
            fout.write(line.replace('\n', ']'))

        for line in fin:
            line = line.rstrip() + ']'
            print >>fout, line
        '''
        for     line in fin:
            fout.write(line.replace(' ', ', ', 1))
