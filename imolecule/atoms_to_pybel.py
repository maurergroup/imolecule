def atoms_to_pybel(aseatoms, infer_bonds=True):
    '''
    Convert ASE Atoms isntance into the json format compatible with

    Args:
        aseatoms : ase.Atoms
            Instance of Atoms from ase package
        infer_bonds : bool
            If `True` bonds will be inferred using openbabel

    Returns:
        mol : dist
            A dictionary with the json format of the molecule
    '''

    from openbabel import pybel
    from openbabel import openbabel as ob
    
    obmol = ob.OBMol()
    obmol.BeginModify()

    for atom in aseatoms:
        obatom = obmol.NewAtom()
        obatom.SetAtomicNum(int(atom.number))
        obatom.SetVector(*atom.position.tolist())

    # If there is no bond data, try to infer them
    if infer_bonds:
        obmol.ConnectTheDots()
        obmol.PerceiveBondOrders()

    # Check for unit cell data
    if any(aseatoms.pbc):
        uc = ob.OBUnitCell()
        uc.SetData(*(ob.vector3(*v) for v in aseatoms.get_cell()))
        uc.SetSpaceGroup('P1')
        obmol.CloneData(uc)
    obmol.EndModify()

    mol = pybel.Molecule(obmol)
    return mol
