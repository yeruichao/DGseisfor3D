import numpy

def vtuxml_head(filename,Nnod,Nele,Ndim,Nv,\
        PointScals=[],PointVects=[],CellScals=[],CellVects=[]):
    NPscal=len(PointScals)
    NPvect=len(PointVects)
    NCscal=len(CellScals)
    NCvect=len(CellVects)
    fh=open(filename+".vtu","wb")
    fh.write('<?xml version="1.0"?>\n')
    fh.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
    fh.write('  <UnstructuredGrid>\n')
    fh.write('    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n' % (Nnod,Nele))
    fh.write('      <Points>\n')
    fh.write('        <DataArray type="Float32" NumberOfComponents="%d" Name="Points" format="appended" offset="0"/>\n' % Ndim)
    fh.write('      </Points>\n')
    fh.write('      <Cells>\n')
    offset=Ndim*Nnod*4+4
    fh.write('        <DataArray type="Int32" Name="connectivity" format="appended" offset="%d"/>\n' % offset)
    offset=offset+Nv*Nele*4+4
    fh.write('        <DataArray type="Int32" Name="offsets" format="appended" offset="%d"/>\n' % offset)
    offset=offset+Nele*4+4
    fh.write('        <DataArray type="Int8" Name="types" format="appended" offset="%d"/>\n' % offset)
    offset=offset+Nele+4
    fh.write('      </Cells>\n')
    if NPscal > 0 or NPvect > 0:
        fh.write('      <PointData>\n')
        for i in range(NPscal):
            fh.write('        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="appended" offset="%d"/>\n' \
                    % (PointScals[i],offset))
            offset=offset+Nnod*4+4
        for i in range(NPvect):
            fh.write('        <DataArray type="Float32" Name="%s" NumberOfComponents="3" format="appended" offset="%d"/>\n' \
                    % (PointVects[i],offset))
            offset=offset+Ndim*Nnode*4+4
        fh.write('      </PointData>\n')
    if NCscal > 0 or NCvect > 0:
        fh.write('      <CellData>\n')
        for i in range(NCscal):
            fh.write('        <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="appended" offset="%d"/>\n' \
                    % (CellScals[i],offset))
            offset=offset+Nele*4+4
        for i in range(NCvect):
            fh.write('        <DataArray type="Float32" Name="%s" NumberOfComponents="3" format="appended" offset="%d"/>\n' \
                    % (CellVects[i],offset))
            offset=offset+Ndim*Nele*4+4
        fh.write('      </CellData>\n')
    fh.write('    </Piece>\n')
    fh.write('  </UnstructuredGrid>\n')
    fh.write('  <AppendedData encoding="raw">\n')
    fh.write(bytearray('_'))
    return fh

def vtuxml_point(fh,nod):
    nod=nod.astype(numpy.float32).T
    fh.write(numpy.int32(nod.size*4))
    for i in range(nod.shape[0]):
        fh.write(numpy.asfortranarray(nod[i]))

def vtuxml_cell(fh,ele):
    ele=ele.astype(numpy.int32).T.copy(order='C')-1
    Nele,Nvtx=ele.shape
    fh.write(numpy.int32(Nvtx*Nele*4))
    for i in range(ele.shape[0]):
        fh.write(ele[i])
    fh.write(numpy.int32(Nele*4))
    l=numpy.arange(1,Nele+1)*Nvtx
    fh.write(l.astype(numpy.int32))
    fh.write(numpy.int32(Nele))
    if Nvtx==4:
        l=numpy.ones(Nele)*10
    elif Nvtx==3:
        l=numpy.ones(Nele)*5
    fh.write(l.astype(numpy.int8))

def vtuxml_var(fh,var):
    var=var.astype(numpy.float32).T
    fh.write(numpy.int32(var.size*4))
    fh.write(numpy.asfortranarray(var))

def vtuxml_end(fh):
    fh.write(chr(10))                                           
    fh.write('  </AppendedData>\n')
    fh.write('</VTKFile>')
    fh.close()

