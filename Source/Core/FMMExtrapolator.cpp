#include "FMMExtrapolator.h"
#include "ScalarField.h"

namespace PhysX {

template <int Dim>
void FMMExtrapolator<Dim>::extrapolate(GridBasedData<Dim, double> &phi, const GridBasedData<Dim, double> &rphi, const GridBasedData<Dim, double> &sphi, const SGridBasedData<Dim, double> &spsi, double theta)
{
	// auto phiView = ScalarFieldView<HermiteIntrpl<Dim, double>>(phi, psi);
	auto phiView = ScalarFieldView<LinearIntrpl<Dim>>(rphi);
	auto sphiView = ScalarFieldView<LinearIntrpl<Dim>>(sphi);
	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		_closestNode[cell] = -1;
		if (sphi[cell] <= 0) 	 _cellType[cell] = 1; 	// solid
		else if (phi[cell] <= 0) _cellType[cell] = 2; 	// liquid
		else if (phi[cell] > 0)  _cellType[cell] = 4; 	// air
	});
	_nodeNormals.clear();
	forEach(_sGrid.nodeGrid, [&](const VectorDi &coord) {
		int nodeidx = _sGrid.nodeGrid.index(coord);
		VectorDd pos = _sGrid.nodeGrid.position(coord);
		int flag = 0;
		// std::cout << "node: (" << coord.transpose() << ") " << '\n';
		for (int i = 0; i < StaggeredGrid<Dim>::numberOfCellEdges(); i++) {
			const int axis = StaggeredGrid<Dim>::cellEdgeAxis(i);
			const VectorDi face = StaggeredGrid<Dim>::nodeFace(coord, i);
			const VectorDi cell0 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 0);
			const VectorDi cell1 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 1);
			if (_sGrid.cellGrid.isValid(cell0) && _sGrid.cellGrid.isValid(cell1)) {
				int faceType = _cellType[cell0] | _cellType[cell1];
				// std::cout << "cell:" << (pos / _sGrid.spacing).transpose() << " (" << cell0.transpose() << ") " << _cellType[cell0] << " (" << cell1.transpose() << ") " << _cellType[cell1] << " " << faceType << '\n';
				if (faceType == 3) 		flag = flag | 1; // solid-liquid
				else if (faceType == 5) flag = flag | 2; // solid-air
				else if (faceType == 6) flag = flag | 4; // air-liquid
			}
		}
		// std::cout << "flag: " << flag << '\n';
		if (flag != 7) return;

		// coord is maybe on the contact line
		RotatedCoord<Dim> rCoord;
		rCoord.dirs[0] = sphiView.gradient(pos);
		rCoord.dirs[1] = phiView.gradient(pos);
		// // check the angle
		// if (rCoord.dirs[0].dot(rCoord.dirs[1]) > 0.9 * rCoord.dirs[0].norm() * rCoord.dirs[1].norm())
		// 	return;
		rCoord.orthogonalize();
		// std::cout << "contact: " << (pos / _sGrid.spacing).transpose() << '\t' << phiView.gradient(pos).transpose() << "\tr:" << rCoord.dirs[1].transpose() << '\n';
		_nodeNormals[nodeidx] = rCoord;
		for (int i = 0; i < StaggeredGrid<Dim>::numberOfCellNodes(); i++) {
			const VectorDi cell = StaggeredGrid<Dim>::nodeCell(coord, i);
			if (_sGrid.cellGrid.isValid(cell) && _mark[cell] && _closestNode[cell] < 0) {
				// std::cout << '\t' << cell.transpose() << '\t' << _cellType[cell];
				_closestNode[cell] = nodeidx;
				_heap.push(HeapElement(0.0, int(_sGrid.cellGrid.index(cell))));
			}
		}
		// std::cout << '\n';
	});

	while (!_heap.empty()) {
		const int cellidx = _heap.top().second;
		const VectorDi cell = _sGrid.cellGrid.coordinate(cellidx);
		const VectorDd pos = _sGrid.cellGrid.position(cell);
		const int nodeidx = _closestNode[cell];
		const VectorDi node = _sGrid.nodeGrid.coordinate(nodeidx);
		const VectorDd nodepos = _sGrid.nodeGrid.position(node);
		const VectorDd disp = pos - nodepos;
		const RotatedCoord<Dim> rCoord = _nodeNormals[nodeidx];
		_virtualPhi[cell] = rCoord.dirs[0].dot(disp) * std::cos(theta) + rCoord.dirs[1].dot(disp) * std::sin(theta);
		// std::cout << "solve: p" << (pos / _sGrid.spacing).transpose() << "\tn" << (nodepos / _sGrid.spacing).transpose() << '\t' << disp.norm() << '\t' << _virtualPhi[cell] << '\t' << _heap.top().first << '\n';
		_heap.pop();
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
			// std::cout << nbCell.transpose() << '\t' << _cellType[nbCell] << '\t' << _mark[nbCell] << '\n';
			if (_sGrid.cellGrid.isValid(nbCell) && _mark[nbCell] && _closestNode[nbCell] < 0) {
				const VectorDd nbpos = _sGrid.cellGrid.position(nbCell);
				const VectorDd nbdisp = pos - nodepos;
				_closestNode[nbCell] = nodeidx;
				_heap.push(HeapElement(nbdisp.norm(), int(_sGrid.cellGrid.index(nbCell))));
			}
		}
	}

	parallelForEach(_sGrid.cellGrid, [&](const VectorDi &cell) {
		if (_mark[cell] && _closestNode[cell] > 0) {
			const RotatedCoord<Dim> rCoord = _nodeNormals[_closestNode[cell]];
			phi[cell] = std::max(_virtualPhi[cell], sphi[cell]);
			// psi[cell] = rCoord.dirs[0] * std::cos(theta) + rCoord.dirs[1] * std::sin(theta);
		}
		if (!_mark[cell] && sphi[cell] < 0) {
			const VectorDd pos = _sGrid.cellGrid.position(cell);
			phi[cell] = -sphi[cell];
			// psi[cell] = -sphiView.gradient(pos);
		}
	});
}

template <int Dim>
void FMMExtrapolator<Dim>::setupMark(const GridBasedData<Dim, double> &phi, GridBasedData<Dim, double> &sphi, int maxSteps)
{
	GridBasedData<Dim, int> _mark2(_mark.grid);
	_mark.setConstant(0);
	_mark2.setConstant(0);
	parallelForEach(phi.grid, [&](const VectorDi &cell) {
		if (sphi[cell] > 0) return;
		if (_mark[cell]) { _mark2[cell] = 1; return; }
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
			if (phi.grid.isValid(nbCell) && sphi[nbCell] > 0 && phi[nbCell] < 0) {
				_mark2[cell] = 1;
				return;
			}
		}
	});
	for (int i = 1; i < maxSteps; ++i) {
		parallelForEach(phi.grid, [&](const VectorDi &cell) {
			if (sphi[cell] > 0) return;
			if (_mark2[cell]) { _mark[cell] = 1; return; }
			for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
				const VectorDi &nbCell = Grid<Dim>::neighbor(cell, i);
				if (phi.grid.isValid(nbCell) && _mark2[nbCell]) {
					_mark[cell] = 1;
					return;
				}
			}
		});
		std::swap(_mark, _mark2);
	}
}

template class FMMExtrapolator<2>;
template class FMMExtrapolator<3>;

}
