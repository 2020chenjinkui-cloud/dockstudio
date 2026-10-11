// @vitest-environment jsdom
import { cleanup, render, screen } from '@testing-library/react';
import { afterEach, describe, expect, it } from 'vitest';
import CrystalRmsdCell from './CrystalRmsdCell';

afterEach(cleanup);

describe('CrystalRmsdCell', () => {
  it('null 显示 —', () => {
    const { container } = render(<CrystalRmsdCell rmsd={null} />);
    expect(container.textContent).toBe('—');
  });

  it('< 2 Å 显示数值与通过徽标', () => {
    render(<CrystalRmsdCell rmsd={1.852} />);
    expect(screen.getByText(/1\.85/)).toBeTruthy();
    expect(screen.getByText('通过')).toBeTruthy();
  });

  it('>= 2 Å 只显示数值无徽标', () => {
    const { container } = render(<CrystalRmsdCell rmsd={3.21} />);
    expect(container.textContent).toContain('3.21');
    expect(screen.queryByText('通过')).toBeNull();
  });
});
