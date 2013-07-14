def iterpose(self, rmsd=0.0001):
        """Iteratively superpose the ensemble until convergence.  Initially, 
        all conformations are aligned with the reference coordinates.  Then 
        mean coordinates are calculated, and are set as the new reference 
        coordinates.  This is repeated until reference coordinates do not 
        change.  This is determined by the value of RMSD between the new and 
        old reference coordinates.  Note that at the end of the iterative 
        procedure the reference coordinate set will be average of conformations
        in the ensemble. 
        
        :arg rmsd: change in reference coordinates to determine convergence,
            default is 0.0001 Å RMSD
        :type rmsd: float"""

        if self._coords is None:
            raise AttributeError('coordinates are not set, use `setCoords`')
        if self._confs is None or len(self._confs) == 0:
            raise AttributeError('conformations are not set, use'
                                    '`addCoordset`')
        LOGGER.info('Starting iterative superposition:')
        LOGGER.timeit('_prody_ensemble')
        rmsdif = 1
        weights = self._weights
        if weights is not None and weights.ndim == 3:
            weightsum = weights.sum(axis=0)
        length = len(self)
        while rmsdif > rmsd:# and step < 1:
                self._superpose()
                newxyz = self._confs.sum(0) / length
                with file("iter_step_%d.coords"%(step), 'w') as outfile:
                        outfile.write("%d %d %d\n"%self._confs.shape)
                        for coordset in self._confs:
                                np.savetxt(outfile,coordset)

                #if weights is None:
                newxyz = self._confs.sum(0) / length
                with file("mean_step_%d.coords"%(step), 'w') as outfile:
                        outfile.write("%d %d %d\n"%(1, newxyz.shape[0], newxyz.shape[1]))
                        np.savetxt(outfile,newxyz)
                #else:
                #newxyz = (self._confs * weights).sum(0) / weightsum
                rmsdif = getRMSD(self._coords, newxyz)
                self._coords = newxyz
                step += 1
                LOGGER.info(('Step #{0}: RMSD difference = '
                      '{1:.4e}').format(step, rmsdif))
        LOGGER.report('Iterative superposition completed in %.2fs.',
                      '_prody_ensemble')
