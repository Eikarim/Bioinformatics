def simpleinteractions(hblist):

	# Output format example:    MET 1N         VAL 17O         2.564

    print ()  #This is the header
    print ("Polar interactions")
    print ('{:13} {:13} {:6} '.format(
            'Atom1','Atom2','Dist (A)')
    )		

    for hb in sorted (hblist,key=lambda i: i[0].get_serial_number()):

        r1 = hb[0].get_parent()
        r2 = hb[1].get_parent()
        
        print ('{:14} {:14} {:6.3f} '.format(
            r1.get_resname()+' '+str(r1.id[1])+hb[0].id,
            r2.get_resname()+' '+str(r2.id[1])+hb[1].id,
            hb[0] - hb[1]
            )
        )

if __name__ == "__main__":
    main()