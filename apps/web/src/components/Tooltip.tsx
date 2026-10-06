interface TooltipProps {
  text: string;
}

export default function Tooltip({ text }: TooltipProps) {
  return (
    <span className="tooltip">
      <span className="tooltip-trigger" tabIndex={0} role="note">
        ?
      </span>
      <span className="tooltip-body">{text}</span>
    </span>
  );
}
